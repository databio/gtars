


/// Attempt to compress counts before writing to bedGraph
pub fn compress_counts(count_results: &mut (Vec<u32>, Vec<i32>), start_position: i32) -> (Vec<u32>, Vec<u32>, Vec<u32>){

    let mut final_starts: Vec<u32> = Vec::new();
    let mut final_ends: Vec<u32> = Vec::new();
    let mut final_counts: Vec<u32> = Vec::new();

    // .0 are the counts, .1 are the positions to track
    let mut previous_count = count_results.0[0];

    let mut previous_start = start_position as u32;
    let mut current_start = previous_start;

    let mut current_end = start_position as u32;


    for (u, i) in count_results.0.iter().zip(count_results.1.iter()) {
        //println!("u: {}, i: {}", u, i);
        let current_count = *u;
        current_end = current_end + 1;

        if current_count != previous_count{
            final_starts.push(current_start);
            final_ends.push(current_end);
            final_counts.push(previous_count);
            current_start = current_end;
            previous_count = current_count;
        } else{
            previous_count = current_count;
        }

    }

    // Must add these lines else we will not get the closing interval (since previous count will be = current count at the close).
    final_starts.push(current_start);
    final_ends.push(current_end);
    final_counts.push(previous_count);

    // println!("Final Starts:{:?}", final_starts);
    // println!("Final Ends:{:?}", final_ends);
    // println!("Final Counts:{:?}", final_counts);

    (final_starts,final_ends, final_counts)


}