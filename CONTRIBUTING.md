## Contributing

### New internal library crate tools

If you'd like to add a new tool, you can do so by creating a new crate at the `gtars` workspace root. This can be done via `cargo new --lib <tool_name>`.

### New public library crate tools

If you want this to be available to users of `gtars`, you can add it to the `gtars` library crate as well. To do so, add your new crate as a dependency in `Cargo.toml`:

```toml
[dependencies]
<tool_name> = { path = "../<tool_name>", optional = true }
```

Then create a new feature to `Cargo.toml`:

```toml
[features]
<tool_name> = ["dep:<tool_name>"]
```

And finally, re-export this tool in `gtars/src/lib.rs`:

```rust
#[cfg(feature = "<tool_name>")]
#[doc(inline)]
pub use <tool_name> as <tool_name>;

```

### New binary crate tools

Finally, if you want to have command-line functionality, you can add it to the `gtars-cli` binary crate. This requires ____ steps:

1. Create a new module inside `gtars-cli/src/main.rs`:

```rust
mod new_tool;
```

2. Use `clap` to define your command-line interface inside `gtars-cli/src/new_tool/cli.rs`:

```rust
use clap::{App, Arg};

pub fn make_new_tool_cli() -> App<'static> {
    App::new("new_tool")
        .about("Does something new")
        .arg(
            Arg::new("input")
                .about("Input file")
                .required(true)
                .index(1),
        )
}
```

3. Write your logic in a wrapper function. This will live inside the `handlers` module of `gtars-cli/src/new_tool`:

```rust
// top of file:
use tool_name::{ ... }

// inside the module:
pub fn new_tool_wrapper() -> Result<(), Box<dyn Error>> {
    // your logic here
}
```

4. Bring this into the `gtars-cli/src/main.rs` file:

```rust
mod new_tool;
```

Please make sure you update the changelog and bump the version number in `Cargo.toml` when you add a new tool.


## Testing

To run the tests, run `cargo test --workspace`.

### Refget tests

The default tests for this module are designed to run quickly on tiny fasta files.
To run the test on a full-scale fasta file, you can look at `test_loading_large_fasta_file`.
This is large test, which is ignored by default, so it doesn't run in the typical `cargo test`. 
To run just this large test on a fasta file, try something like this:

```
FASTA_PATH=tests/data/subset.fa.gz cargo test tests::test_loading_large_fasta_file -- --nocapture --ignored
FASTA_PATH=`refgenie seek test/fasta` cargo test tests::test_loading_large_fasta_file -- --nocapture --ignored
```
