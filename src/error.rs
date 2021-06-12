use std::fmt;

#[derive(Debug)]
pub struct CastError<'a, F> {
    from: &'a F,
    to: String,
}

impl<'a, F> CastError<'a, F> {
    pub fn new(from: &'a F, to: &str) -> Self {
        Self {
            from,
            to: to.into(),
        }
    }
}

impl<'a, F> fmt::Display for CastError<'a, F>
where
    F: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[ERROR] Could not cast a {:?} to a {:?}",
            self.from, self.to
        )
    }
}

impl<'a, F: fmt::Debug> std::error::Error for CastError<'a, F> {}
