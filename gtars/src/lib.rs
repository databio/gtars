#[cfg(feature = "core")]
#[doc(inline)]
pub use gtars_core as core;

#[cfg(feature = "tokenizers")]
#[doc(inline)]
pub use tokenizers as tokenizers;

#[cfg(feature = "io")]
#[doc(inline)]
pub use io as io;

#[cfg(feature = "refget")]
#[doc(inline)]
pub use refget as refget;

#[cfg(feature = "overlaprs")]
#[doc(inline)]
pub use overlaprs as overlaprs;

#[cfg(feature = "uniwig")]
#[doc(inline)]
pub use uniwig as uniwig;

#[cfg(feature = "igd")]
#[doc(inline)]
pub use igd as igd;

#[cfg(feature = "bbcache")]
#[doc(inline)]
pub use bbcache as bbcache;

#[cfg(feature = "scoring")]
#[doc(inline)]
pub use scoring as scoring;

#[cfg(feature = "fragsplit")]
#[doc(inline)]
pub use fragsplit as fragsplit;
