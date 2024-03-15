use clap::ArgMatches;

pub mod cli;


pub fn read_bed_map(){

}

pub fn read_bed_bec(){

}

pub fn parse_bed_file(line: &str) -> Option<(String, i32, i32, Option<String>)> {
    // TODO Eventually refactor all bed file parsing to a single shared function

    let mut iter = line.split('\t');
    let mut ctg: Option<String> = None;
    let mut st = -1;
    let mut en = -1;
    let mut r: Option<String> = None;
    let mut i = 0;

    let iter = iter.by_ref();

    while let Some(item) = iter.next()  {
        match i {
            0 => ctg = Some(item.to_owned()),
            1 => st = item.parse().unwrap_or(-1),
            2 => {
                en = item.parse().unwrap_or(-1);
                r = iter.next().map(|x| x.to_owned());
            }
            _ => break,
        }
        i += 1;
        if i == 3 || iter.next().is_none() {
            break;
        }
    }

    Some((ctg?, st, en, r))

}


pub fn run_uniwig(matches: &ArgMatches) {
    println!("Im running. Here are the arguments: {:?}", matches)
}

pub mod consts {
    pub const UNIWIG_CMD: &str = "uniwig";

}