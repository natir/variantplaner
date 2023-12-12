//! Function required to compute variant id

/* std use */

/* crate use */

/* polars use */
use polars_core::prelude::*;
use pyo3_polars::derive::polars_expr;

#[inline(always)]
pub(crate) fn nuc2bit(nuc: u8) -> u64 {
    (nuc as u64 >> 1) & 0b11
}

#[inline(always)]
pub(crate) fn seq2bit(seq: &[u8]) -> u64 {
    let mut two_bit = 0;

    for nuc in seq {
        two_bit <<= 2;
        two_bit |= nuc2bit(*nuc)
    }

    two_bit
}

fn local_compute(
    real_pos: &UInt64Chunked,
    ref_seq: &Utf8Chunked,
    alt_seq: &Utf8Chunked,
    max_pos: u64,
) -> PolarsResult<Series> {
    let pos_mov = max_pos.leading_zeros() as u64 - 1;
    let sep_len = (pos_mov as f64 / 2.0).floor().log2().ceil() as u64;
    let sep_mov = pos_mov - sep_len;
    let nuc_len_max = sep_mov / 2;
    let hasher = ahash::RandomState::with_seeds(42, 42, 42, 42);
    let mut key = Vec::with_capacity(128);

    let out: ChunkedArray<UInt64Type> = real_pos
        .into_iter()
        .zip(ref_seq)
        .zip(alt_seq)
        .map(|((p, r), a)| match (p, r, a) {
            (Some(p), Some(r), Some(a)) => {
                let mut hash = 0;
                if r.len() + a.len() > nuc_len_max as usize {
                    key.clear();

                    key.extend(p.to_be_bytes());
                    key.extend(r.as_bytes());
                    key.extend(a.as_bytes());

                    hash = (1 << 63) | (hasher.hash_one(&key) >> 1);
                } else {
                    hash |= p << pos_mov;
                    hash |= (r.len() as u64) << sep_mov;
                    hash |= seq2bit(r.as_bytes()) << (a.len() * 2);
                    hash |= seq2bit(a.as_bytes());
                }
                Some(hash)
            }
            _ => None,
        })
        .collect();

    Ok(out.into_series())
}

#[polars_expr(output_type=UInt64)]
fn compute(inputs: &Vec<Series>) -> PolarsResult<Series> {
    let real_pos = inputs[0].u64()?;
    let ref_seq = inputs[1].utf8()?;
    let alt_seq = inputs[2].utf8()?;
    let max_pos = inputs[3].cast(&DataType::UInt64)?.u64()?.get(0).unwrap();

    local_compute(real_pos, ref_seq, alt_seq, max_pos)
}

fn local_part(id: &UInt64Chunked) -> PolarsResult<Series> {
    Ok(id
        .into_iter()
        .map(|i| match i {
            Some(i) if i >> 63 == 0b1 => Some(255),
            Some(i) => Some((i << 1) >> 56),
            _ => None,
        })
        .collect())
}

#[polars_expr(output_type=UInt8)]
fn partition(inputs: &Vec<Series>) -> PolarsResult<Series> {
    let id = inputs[0].u64()?;

    local_part(id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nuc2bit_() {
        assert_eq!(nuc2bit(b'A'), 0b00);
        assert_eq!(nuc2bit(b'C'), 0b01);
        assert_eq!(nuc2bit(b'T'), 0b10);
        assert_eq!(nuc2bit(b'G'), 0b11);
        assert_eq!(nuc2bit(b'N'), 0b11);

        assert_eq!(nuc2bit(b'a'), 0b00);
        assert_eq!(nuc2bit(b'c'), 0b01);
        assert_eq!(nuc2bit(b't'), 0b10);
        assert_eq!(nuc2bit(b'g'), 0b11);
        assert_eq!(nuc2bit(b'n'), 0b11);
    }

    #[test]
    fn seq2bit_() {
        assert_eq!(seq2bit(b"ACTGN"), 0b0001101111);
    }

    #[test]
    fn compute_id() {
        let mut real_pos = UInt64Chunked::new_vec(
            "real_pos",
            vec![10, 50, 110, 326512443305, 326512443305, 224],
        );
        let mut ref_seq = Utf8Chunked::new("ref", vec!["A", "C", "T", "G", "GA", "CATGAGCGGACTG"]);
        let mut alt_seq = Utf8Chunked::new("alt", vec!["G", "T", "C", "A", "", "AC"]);

        real_pos.extend(&UInt64Chunked::full_null("", 1));
        ref_seq.extend(&Utf8Chunked::full_null("", 1));
        alt_seq.extend(&Utf8Chunked::full_null("", 1));

        let id = local_compute(&real_pos, &ref_seq, &alt_seq, 326512443305).unwrap();

        assert_eq!(
            id,
            Series::from_any_values_and_dtype(
                "",
                &[
                    AnyValue::UInt64(168820739),
                    AnyValue::UInt64(839909382),
                    AnyValue::UInt64(1846542345),
                    AnyValue::UInt64(5477969788016787468),
                    AnyValue::UInt64(5477969788017836044),
                    AnyValue::UInt64(13756669010756803764),
                    AnyValue::Null,
                ],
                &DataType::UInt64,
                false
            )
            .unwrap()
        );
    }

    #[test]
    fn compute_part() {
        let mut real_pos = UInt64Chunked::new_vec(
            "real_pos",
            vec![10, 50, 110, 326512443305, 326512443305, 224],
        );
        let mut ref_seq = Utf8Chunked::new("ref", vec!["A", "C", "T", "G", "GA", "CATGAGCGGACTG"]);
        let mut alt_seq = Utf8Chunked::new("alt", vec!["G", "T", "C", "A", "", "AC"]);

        real_pos.extend(&UInt64Chunked::full_null("", 1));
        ref_seq.extend(&Utf8Chunked::full_null("", 1));
        alt_seq.extend(&Utf8Chunked::full_null("", 1));

        let id = local_compute(&real_pos, &ref_seq, &alt_seq, 326512443305).unwrap();
        let partition = local_part(id.u64().unwrap()).unwrap();

        assert_eq!(
            partition,
            Series::from_any_values_and_dtype(
                "",
                &[
                    AnyValue::UInt8(0),
                    AnyValue::UInt8(0),
                    AnyValue::UInt8(0),
                    AnyValue::UInt8(152),
                    AnyValue::UInt8(152),
                    AnyValue::UInt8(255),
                    AnyValue::Null,
                ],
                &DataType::UInt8,
                false
            )
            .unwrap()
        );
    }
}
