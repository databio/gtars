//! `rgblast search`: rgblast's own flags are not forwarded to the engine, wherever they sit.

mod common;

use std::process::Command;

#[test]
fn quiet_after_blastn_flags_stays_out_of_the_engine_args() {
    let t = common::tempdir();
    let fake = t.path().join("fake_blastn.sh");
    let seen = t.path().join("seen");
    std::fs::write(
        &fake,
        format!(
            "#!/bin/sh\ncase \"$1\" in\n  -version) echo 'blastn: 2.17.0+'; exit 0 ;;\nesac\n\
             echo \"$@\" > {}\nfor a in \"$@\"; do case \"$a\" in --quiet|--no-extra-fields|--meta) exit 2 ;; esac; done\n",
            seen.display()
        ),
    )
    .unwrap();
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    let o = Command::new(env!("CARGO_BIN_EXE_rgblast"))
        .args(["search", "--engine", "blastn", "--db"])
        .arg(t.path().join("db"))
        .arg("--blastn")
        .arg(&fake)
        .args([
            "-query",
            "q.fa",
            "-outfmt",
            "6",
            "--quiet",
            "--no-extra-fields",
        ])
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert_eq!(o.status.code(), Some(0), "{err}");
    let got = std::fs::read_to_string(&seen).unwrap();
    assert!(got.contains("-query q.fa -outfmt 6"), "{got}");
    assert!(!got.contains("--quiet"), "{got}");
}
