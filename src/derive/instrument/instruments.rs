//! Knowledge about which instrument ids map to which machine types.

use std::collections::{HashMap, HashSet};

/// Encapsulates the knowledge we currently have on which instrument name patterns map
/// to which machine types as a [`HashMap`].
pub fn build_instrument_lookup_table() -> HashMap<&'static str, HashSet<&'static str>> {
    HashMap::from([
        ("^HWI-M[0-9]{4}$", HashSet::from(["MiSeq"])),
        ("^HWUSI", HashSet::from(["Genome Analyzer IIx"])),
        ("^M[0-9]{5}$", HashSet::from(["MiSeq"])),
        ("^HWI-C[0-9]{5}$", HashSet::from(["HiSeq 1500"])),
        ("^C[0-9]{5}$", HashSet::from(["HiSeq 1500"])),
        (
            "^HWI-ST[0-9]{3,5}(_[0-9]{9})?$",
            HashSet::from(["HiSeq 2000"]),
        ),
        (
            "^HWI-D[0-9]{5}$",
            HashSet::from(["HiSeq 2000", "HiSeq 2500"]),
        ),
        ("^A[0-9]{5}$", HashSet::from(["NovaSeq"])),
        ("^D[0-9]{5}$", HashSet::from(["HiSeq 2500"])),
        ("^J[0-9]{5}$", HashSet::from(["HiSeq 3000"])),
        ("^K[0-9]{5}$", HashSet::from(["HiSeq 3000", "HiSeq 4000"])),
        ("^E[0-9]{5}$", HashSet::from(["HiSeq X"])),
        ("^N[0-9]{5}$", HashSet::from(["NextSeq"])),
        ("^NB[0-9]{6}$", HashSet::from(["NextSeq"])),
        ("^NS[0-9]{6}$", HashSet::from(["NextSeq"])),
        ("^MN[0-9]{5}$", HashSet::from(["MiniSeq"])),
    ])
}
