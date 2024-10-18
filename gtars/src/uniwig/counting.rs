/// This function is a more direct port of smoothFixedStartEndBW from uniwig written in CPP.
/// It allows the user to accumulate reads of either starts or ends.
/// Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
/// the level of smoothing.
/// counts are reported over a stepsize (with a default of stepsize = 1).
/// Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.
#[allow(unused_variables)]
pub fn start_end_counts(
    starts_vector: &Vec<(i32, i32)>,
    chrom_size: i32,
    smoothsize: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {
    //let vin_iter = starts_vector.iter();

    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count: i32 = 0;

    let mut coordinate_value: (i32, i32);
    let mut prev_coordinate_value = 0;

    let mut adjusted_start_site: (i32, i32);
    let mut current_end_site: (i32, i32);

    let mut collected_end_sites: Vec<(i32, i32)> = Vec::new();

    adjusted_start_site = starts_vector[0].clone(); // get first coordinate position

    adjusted_start_site.0 = adjusted_start_site.0 - smoothsize;

    current_end_site = adjusted_start_site;
    current_end_site.0 = adjusted_start_site.0 + 1 + smoothsize * 2;

    if adjusted_start_site.0 < 1 {
        adjusted_start_site.0 = 1;
    }

    while coordinate_position < adjusted_start_site.0 {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    for (index, coord) in starts_vector.iter().enumerate().skip(0) {
        coordinate_value = *coord;

        adjusted_start_site = coordinate_value;
        adjusted_start_site.0 = coordinate_value.0 - smoothsize;

        let current_score = adjusted_start_site.1;

        count += current_score;

        if adjusted_start_site.0 < 1 {
            adjusted_start_site.0 = 1;
        }

        let current_index = index;

        let mut new_end_site = adjusted_start_site;
        new_end_site.0 = adjusted_start_site.0 + 1 + smoothsize * 2;
        collected_end_sites.push(new_end_site);

        if adjusted_start_site.0 == prev_coordinate_value {
            continue;
        }

        while coordinate_position < adjusted_start_site.0 {
            while current_end_site.0 == coordinate_position {
                count = count - current_score;

                if collected_end_sites.last() == None {
                    current_end_site.0 = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count as u32);
                v_coordinate_positions.push(coordinate_position);
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = adjusted_start_site.0;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.
                       // this is because the code above subtracts twice during the INITIAL end site closure. So we are missing one count and need to make it up else we go negative.

    while coordinate_position < chrom_size {
        // Apply a bound to push the final coordinates otherwise it will become truncated.

        while current_end_site.0 == coordinate_position {
            let current_score = adjusted_start_site.1;
            count = count - current_score;

            if collected_end_sites.last() == None {
                current_end_site.0 = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count as u32);
            v_coordinate_positions.push(coordinate_position);
        }

        coordinate_position = coordinate_position + 1;
    }

    (v_coord_counts, v_coordinate_positions)
}

/// This function is a more direct port of fixedCoreBW from uniwig written in CPP
/// It allows the user to accumulate reads across paired starts and ends.
/// Counts occur between a start coordinate (cutSite) and an end site (endSite) where the endsite is determined based on
/// the paired ends.
/// Counts are reported over a stepsize (with a default of stepsize = 1)
/// Unlike the original function, it does not write to disk in chunks. it simply returns a vector of accumulated reads.
#[allow(unused_variables)]
pub fn core_counts(
    starts_vector: &Vec<(i32, i32)>,
    ends_vector: &Vec<(i32, i32)>,
    chrom_size: i32,
    stepsize: i32,
) -> (Vec<u32>, Vec<i32>) {
    let mut v_coordinate_positions: Vec<i32> = Vec::new(); // these are the final coordinates after any adjustments
    let mut v_coord_counts: Vec<u32> = Vec::new(); // u8 stores 0:255 This may be insufficient. u16 max is 65535

    let mut coordinate_position = 1;

    let mut count = 0;

    let mut coordinate_value: (i32, i32);
    let mut prev_coordinate_value = 0;

    let mut current_start_site: (i32, i32);
    let mut current_end_site: (i32, i32);

    let mut collected_end_sites: Vec<(i32, i32)> = Vec::new();

    current_start_site = starts_vector[0].clone(); // get first coordinate position
    current_end_site = ends_vector[0];

    if current_start_site.0 < 1 {
        current_start_site.0 = 1;
    }

    while coordinate_position < current_start_site.0 {
        // Just skip until we reach the initial adjusted start position
        // Note that this function will not return 0s at locations before the initial start site
        coordinate_position = coordinate_position + stepsize;
    }

    for (index, coord) in starts_vector.iter().enumerate().skip(0) {
        coordinate_value = *coord;

        current_start_site = coordinate_value;

        let current_score = current_start_site.1;
        count += current_score;

        if current_start_site.0 < 1 {
            current_start_site.0 = 1;
        }

        let current_index = index;

        collected_end_sites.push(ends_vector[current_index]);

        if current_start_site.0 == prev_coordinate_value {
            continue;
        }

        while coordinate_position < current_start_site.0 {
            while current_end_site.0 == coordinate_position {
                count = count - current_score;

                if collected_end_sites.last() == None {
                    current_end_site.0 = 0;
                } else {
                    current_end_site = collected_end_sites.remove(0)
                }
            }

            if coordinate_position % stepsize == 0 {
                // Step size defaults to 1, so report every value
                v_coord_counts.push(count as u32);
                v_coordinate_positions.push(coordinate_position);
            }

            coordinate_position = coordinate_position + 1;
        }

        prev_coordinate_value = current_start_site.0;
    }

    count = count + 1; // We must add 1 extra value here so that our calculation during the tail as we close out the end sites does not go negative.

    while coordinate_position < chrom_size {
        while current_end_site.0 == coordinate_position {
            let current_score = current_start_site.1;
            count = count - current_score;

            if collected_end_sites.last() == None {
                current_end_site.0 = 0;
            } else {
                current_end_site = collected_end_sites.remove(0)
            }
        }

        if coordinate_position % stepsize == 0 {
            // Step size defaults to 1, so report every value
            v_coord_counts.push(count as u32);
            v_coordinate_positions.push(coordinate_position);
        }

        coordinate_position = coordinate_position + 1;
    }

    (v_coord_counts, v_coordinate_positions)
}
