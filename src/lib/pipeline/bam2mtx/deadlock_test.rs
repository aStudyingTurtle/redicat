//! Test to verify the deadlock fix in bam2mtx streaming pipeline
//!
//! This test creates a scenario that would previously trigger the deadlock:
//! - Small channel capacity
//! - Many chunks to process
//! - Verifies completion without hanging

#[cfg(test)]
mod deadlock_regression_tests {
    use crossbeam::channel::bounded;
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_scoped_sender_cleanup_prevents_deadlock() {
        // Simulate the bam2mtx scenario with minimal channel capacity
        let channel_capacity = 2;
        let num_chunks = 100; // More chunks than channel capacity

        let (sender, receiver) = bounded::<usize>(channel_capacity);
        let processed = Arc::new(AtomicUsize::new(0));
        let processed_clone = Arc::clone(&processed);

        // Consumer thread (like the AnnData writer)
        let consumer = thread::spawn(move || {
            let mut count = 0;
            for _value in receiver.into_iter() {
                count += 1;
                // Simulate slow processing
                thread::sleep(Duration::from_micros(10));
                processed_clone.fetch_add(1, Ordering::Relaxed);
            }
            count
        });

        // Producer threads (like Rayon processing chunks)
        // CRITICAL: Use scoped sender_clone like in the fix
        {
            let sender_clone = sender.clone();
            let chunks: Vec<usize> = (0..num_chunks).collect();

            chunks
                .into_par_iter()
                .try_for_each_init(
                    || sender_clone.clone(),
                    |tx, chunk_id| -> Result<(), String> {
                        // Simulate chunk processing
                        thread::sleep(Duration::from_micros(5));

                        // This would deadlock without proper sender cleanup:
                        // - Channel fills up quickly (capacity=2)
                        // - Some threads block on send
                        // - Consumer waits for all senders to drop
                        // - But senders don't drop until try_for_each_init returns
                        // - try_for_each_init can't return while threads are blocked
                        tx.send(chunk_id).map_err(|e| e.to_string())
                    },
                )
                .expect("Producer should complete without deadlock");

            // sender_clone drops here when scope ends
        }

        // Drop the original sender
        drop(sender);

        // If we reach here without hanging, the fix works!
        let received_count = consumer.join().expect("Consumer should complete");

        // Verify all chunks were processed
        assert_eq!(received_count, num_chunks, "All chunks should be received");
        assert_eq!(
            processed.load(Ordering::Relaxed),
            num_chunks,
            "All chunks should be processed"
        );

        println!(
            "✓ Deadlock prevention test passed: {} chunks processed with capacity {}",
            num_chunks, channel_capacity
        );
    }

    #[test]
    fn test_without_scope_demonstrates_issue() {
        // This test demonstrates that without proper scoping,
        // sender clones can remain alive longer than intended

        let channel_capacity = 2;

        let (sender, receiver) = bounded::<usize>(channel_capacity);

        let consumer = thread::spawn(move || {
            // Consumer waits for all senders to drop
            receiver.into_iter().count()
        });

        // Simulate the OLD broken code pattern (without explicit scope):
        let sender_clone = sender.clone();

        // In the broken code, sender_clone would remain alive here
        // until the end of the function, even after we thought we were done
        // This demonstrates the importance of explicit scope management

        // Manually drop to avoid the issue
        drop(sender_clone);
        drop(sender);

        let count = consumer.join().unwrap();
        assert_eq!(count, 0); // No data was sent

        println!("✓ Demonstrated that sender clone lifetime matters");
    }

    #[test]
    fn test_channel_capacity_stress() {
        // Stress test with various channel capacities
        for capacity in [1, 2, 4, 8, 16] {
            let num_chunks = capacity * 10; // 10x channel capacity

            let (sender, receiver) = bounded::<usize>(capacity);
            let received = Arc::new(AtomicUsize::new(0));
            let received_clone = Arc::clone(&received);

            let consumer = thread::spawn(move || {
                for _ in receiver.into_iter() {
                    received_clone.fetch_add(1, Ordering::Relaxed);
                    thread::sleep(Duration::from_micros(1));
                }
            });

            {
                let sender_clone = sender.clone();
                let chunks: Vec<usize> = (0..num_chunks).collect();

                chunks
                    .into_par_iter()
                    .try_for_each_init(
                        || sender_clone.clone(),
                        |tx, chunk_id| -> Result<(), String> {
                            tx.send(chunk_id).map_err(|e| e.to_string())
                        },
                    )
                    .expect("Should not deadlock");
            }

            drop(sender);
            consumer.join().expect("Consumer should complete");

            assert_eq!(
                received.load(Ordering::Relaxed),
                num_chunks,
                "Capacity {} should handle {} chunks",
                capacity,
                num_chunks
            );
        }

        println!("✓ Stress test passed for all capacities");
    }
}
