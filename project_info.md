# REDICAT 项目概览

本文件总结 `src/` 目录下的核心代码结构、关键功能模块以及它们之间的协作方式，帮助快速理解 REDICAT（RNA Editing Cellular Assessment Toolkit）的整体实现。

## 顶层结构

- `main.rs`：命令行入口，基于 `structopt` 解析子命令 `bulk`、`bam2mtx`、`call`，并启用 `env_logger`。旧 `preprocess` 子命令已移除，其核心过滤逻辑被整合进 `bulk` 与 `bam2mtx`。
- `commands/`：定义各个 CLI 子命令的参数与执行逻辑。
  - `base_depth.rs` (`bulk`)：对 BAM/CRAM 进行 pileup，输出每个位点的碱基深度、插入/缺失统计等。利用 `redicat_lib::par_granges` 并行，并在命令层统一处理 MAPQ、条形码、UB 过滤。
  - `bam2mtx.rs`：将单细胞 BAM 转换为稀疏矩阵（AnnData/H5AD），核心流程包括载入条形码、读取/生成位点 TSV、分块并行处理 pileup、UMI 去重及 AnnData 输出；新增 `--two-pass` 参数自动调用 `bulk` 获取第一轮位点集合。
  - `call.rs`：RNA 编辑检测主流程。依次执行输入校验、AnnData 读取、REDIportal 注释、参考/突变矩阵构建、Cell Editing Index 计算以及 mismatch 统计。
- `lib/`：共享库代码，供 CLI 与其他组件复用。
  - `bam2mtx/`：单细胞矩阵构建的底层实现。
    - `anndata_output.rs`：高性能稀疏矩阵构建与 AnnData 写入，使用 `nalgebra_sparse`、`rayon`、`smallvec` 等进行并行化和内存优化。
    - `barcode.rs`：条形码白名单加载与校验。
    - `processor.rs`：核心 BAM pileup 处理，完成 UMI 共识、按细胞/链向累积碱基计数。
    - `region_processor.rs`：把 BAM 处理包装成 `par_granges::RegionProcessor` 以按区域并行。
    - `utils.rs`：位点 TSV 读取、分块、进度条等辅助函数。
  - `call/`：RNA 编辑分析管线组件。
    - `anndata_ops.rs`：AnnData 读写、矩阵层合并、元数据维护。
    - `base_matrix.rs`、`editing_analysis.rs`、`editing.rs`、`reference_genome.rs`、`sparse_ops.rs`、`validation.rs`：分别负责基线矩阵操作、编辑事件分析、编辑类型定义、参考基因组访问、稀疏运算及输入参数校验。
  - 其他模块：
    - `par_granges.rs`：高吞吐并行区域遍历框架。
    - `position/`：pileup 位点/区间的数据结构。
    - `read_filter.rs`：读过滤接口与默认实现（基于 MAPQ 与 SAM flag）。
    - `utils.rs`：通用工具（线程/文件/压缩处理等）。

## 关键流程与协作

### Bulk（Base Depth）
1. 使用 `ParGranges` 把基因组分块并行处理。
2. `BaseProcessor` 在每个区域执行 pileup -> `PileupPosition`，通过自定义 `BulkReadFilter` 同时应用 MAPQ、CB 白名单、UB 合法性与 base quality 过滤。
3. 输出 TSV/GZIP，供后续 `bam2mtx` 或分析使用；`bam2mtx --two-pass` 会自动触发该流程生成 `1pass.tsv.gz`。

### BAM → Matrix
1. 解析 CLI 参数，加载条形码白名单（`BarcodeProcessor`）。
2. 若启用 `--two-pass` 自动执行 `bulk` 第一轮生成位点列表；随后读取 TSV，并按染色体/位置排序与分块。
3. `OptimizedChunkProcessor` 针对每个分块重新打开索引 BAM，遍历 pileup，过滤 reads（MAPQ/baseQ/CB/UMI），并执行 UMI 共识。
4. 合并结果后交给 `AnnDataConverter` 构建稀疏矩阵：
   - 收集唯一细胞/位点标识。
   - 创建索引映射，生成前向/反向链四个碱基矩阵（共 8 层）。
   - 输出 AnnData/H5AD。

### Call（RNA 编辑检测）
1. 校验输入文件与配置。
2. 读取单细胞 AnnData，加载参考基因组与 REDIportal 编辑位点。
3. 注释并筛选候选位点，构建参考/突变矩阵。
4. 计算 Cell Editing Index 与 mismatch 统计，写回 AnnData。

## 并行与性能优化
- 全局使用 `rayon` 提供多线程支持（按需自适应线程池）。
- `ParGranges` 通过 channel + super-chunk 管线实现 IO 与计算流水线。
- `bam2mtx` 在构建稀疏矩阵时采用线程局部缓冲（`UnsafeCell` + smallvec），最终合并为 CSR。
- `AnnData` 输出支持压缩与批处理，减少内存峰值。

## 依赖亮点
- `rust-htslib`：BAM/CRAM 读取、pileup、索引处理。
- `nalgebra_sparse`：稀疏矩阵表示与运算。
- `anndata`/`anndata_hdf5`：H5AD 文件存取。
- `rayon`、`crossbeam`：并行任务调度。
- `flate2`、`gzp`：压缩 IO。
- `polars`：AnnData 元数据的 DataFrame 处理。

## 测试与示例
- 部分单元测试覆盖 CLI 参数解析与基本处理路径（例如 `bam2mtx::tests::test_bam2mtx_args`）。
- `docs/examples.md` 提供典型命令行调用示例。
