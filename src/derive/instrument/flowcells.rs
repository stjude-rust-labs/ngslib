//! Knowledge about which flowcells map to which machine types.

use std::collections::{HashMap, HashSet};

/// Encapsulates the knowledge we currently have on which flowcell patterns map
/// to which machine types as a [`HashMap`].
pub fn build_flowcell_lookup_table() -> HashMap<&'static str, HashSet<&'static str>> {
    HashMap::from([
        (
            // High Output (8-lane) v4 flow cell
            "^C[A-Z0-9]{4}ANXX$",
            HashSet::from(["HiSeq 1500", "HiSeq 2000", "HiSeq 2500"]),
        ),
        (
            // High Output (8-lane) v3 flow cell
            "^C[A-Z0-9]{4}ACXX$",
            HashSet::from(["HiSeq 1000", "HiSeq 1500", "HiSeq 2000", "HiSeq 2500"]),
        ),
        (
            // High Output (8-lane) v3 flow cell
            "^D[A-Z0-9]{4}ACXX$",
            HashSet::from(["HiSeq 1000", "HiSeq 1500", "HiSeq 2000", "HiSeq 2500"]),
        ),
        (
            // Rapid Run (2-lane) v1 flow cell
            "^H[A-Z0-9]{4}ADXX$",
            HashSet::from(["HiSeq 1500", "HiSeq 2000", "HiSeq 2500"]),
        ),
        (
            // Rapid Run (2-lane) v2 flow cell
            "^H[A-Z0-9]{4}BCXX$",
            HashSet::from(["HiSeq 1500", "HiSeq 2500"]),
        ),
        (
            // Rapid Run (2-lane) v2 flow cell
            "^H[A-Z0-9]{4}BCXY$",
            HashSet::from(["HiSeq 1500", "HiSeq 2500"]),
        ),
        (
            // (8-lane) v1 flow cell
            "^H[A-Z0-9]{4}BBXX$",
            HashSet::from(["HiSeq 4000"]),
        ),
        (
            // (8-lane) v1 flow cell
            "^H[A-Z0-9]{4}BBXY$",
            HashSet::from(["HiSeq 4000"]),
        ),
        (
            // (8-lane) flow cell
            "^H[A-Z0-9]{4}CCXX$",
            HashSet::from(["HiSeq X"]),
        ),
        (
            // (8-lane) flow cell
            "^H[A-Z0-9]{4}CCXY$",
            HashSet::from(["HiSeq X"]),
        ),
        (
            // (8-lane) flow cell
            "^H[A-Z0-9]{4}ALXX$",
            HashSet::from(["HiSeq X"]),
        ),
        (
            // High output flow cell
            "^H[A-Z0-9]{4}BGX[A-Z,0-9]$",
            HashSet::from(["NextSeq"]),
        ),
        (
            // Mid output flow cell
            "^H[A-Z0-9]{4}AFXX$",
            HashSet::from(["NextSeq"]),
        ),
        (
            // S1 flow cell
            "^H[A-Z0-9]{5}RXX$",
            HashSet::from(["NovaSeq"]),
        ),
        (
            // SP flow cell
            "^H[A-Z0-9]{5}RXX$",
            HashSet::from(["NovaSeq"]),
        ),
        (
            // S2 flow cell
            "^H[A-Z0-9]{5}MXX$",
            HashSet::from(["NovaSeq"]),
        ),
        (
            // S4 flow cell
            "^H[A-Z0-9]{5}SXX$",
            HashSet::from(["NovaSeq"]),
        ),
        (
            // MiSeq flow cell
            "^A[A-Z0-9]{4}$",
            HashSet::from(["MiSeq"]),
        ),
        (
            // MiSeq flow cell
            "^B[A-Z0-9]{4}$",
            HashSet::from(["MiSeq"]),
        ),
        (
            // MiSeq nano flow cell
            "^D[A-Z0-9]{4}$",
            HashSet::from(["MiSeq"]),
        ),
        // (
        //     // Unknown HiSeq flow cell from SJ data
        //     "^D[A-Z0-9]{4}$",
        //     HashSet::from(["HiSeq 2000", "HiSeq 2500"]),
        // ),
        (
            // MiSeq micro flow cell
            "^G[A-Z0-9]{4}$",
            HashSet::from(["MiSeq"]),
        ),
    ])
}
