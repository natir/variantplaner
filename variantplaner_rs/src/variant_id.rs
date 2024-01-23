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

#[inline(always)]
pub(crate) fn ref_alt_space_usage(refs: &[u8], alts: &[u8]) -> u64 {
    (refs.len() as f64).log2().ceil() as u64 + (alts.len() as u64 * 2)
}

fn local_compute(
    real_pos: &UInt64Chunked,
    ref_seq: &Utf8Chunked,
    alt_seq: &Utf8Chunked,
    max_pos: u64,
) -> PolarsResult<Series> {
    let pos_mov = max_pos.leading_zeros() as u64 - 1;
    let hasher = ahash::RandomState::with_seeds(42, 42, 42, 42);
    let mut key = Vec::with_capacity(128);

    println!("max_pos, pos_mov {}", pos_mov);

    let out: ChunkedArray<UInt64Type> = real_pos
        .into_iter()
        .zip(ref_seq)
        .zip(alt_seq)
        .map(|((p, r), a)| match (p, r, a) {
            (Some(p), Some(r), Some(a)) => {
                let mut hash = 0;
                if ref_alt_space_usage(r.as_bytes(), a.as_bytes()) > pos_mov {
                    key.clear();

                    key.extend(p.to_be_bytes());
                    key.extend(r.as_bytes());
                    key.extend(a.as_bytes());

                    hash = (1 << 63) | (hasher.hash_one(&key) >> 1);
                } else {
                    hash |= p << pos_mov;
                    hash |= (r.len() as u64) << (a.len() * 2);
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
        let mut ref_seq = Utf8Chunked::new("ref", vec!["A", "C", "T", "G", "GA", "AC"]);
        let mut alt_seq = Utf8Chunked::new("alt", vec!["G", "T", "C", "A", "", "CATGAGCGGACTG"]);

        real_pos.extend(&UInt64Chunked::full_null("", 1));
        ref_seq.extend(&Utf8Chunked::full_null("", 1));
        alt_seq.extend(&Utf8Chunked::full_null("", 1));

        let id = local_compute(&real_pos, &ref_seq, &alt_seq, 326512443305).unwrap();

        assert_eq!(
            id,
            Series::from_any_values_and_dtype(
                "",
                &[
                    AnyValue::UInt64(167772167),
                    AnyValue::UInt64(838860806),
                    AnyValue::UInt64(1845493765),
                    AnyValue::UInt64(5477969788015738884),
                    AnyValue::UInt64(5477969788015738882),
                    AnyValue::UInt64(17294972461992989018),
                    AnyValue::Null,
                ],
                &DataType::UInt64,
                false
            )
            .unwrap()
        );
    }

    #[test]
    fn homopolymer_colision() {
        let real_pos = UInt64Chunked::new_vec("real_pos", vec![110; 400]);

        let mut refs = vec!["A"; 4];
        refs.extend(vec!["AA"; 20]);
        let ref_seq = Utf8Chunked::new("ref", refs.clone());

        let alts = vec![
            "A", "C", "G", "T", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG",
            "GT", "TA", "TC", "TG", "TT", "A", "C", "G", "T",
        ];

        let alt_seq = Utf8Chunked::new("alt", alts.clone());

        let mut ids: Vec<u64> = local_compute(&real_pos, &ref_seq, &alt_seq, 326512443305)
            .unwrap()
            .u64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        ids.sort_unstable();
        ids.dedup();

        assert_eq!(ids.len(), 24);
    }

    #[test]
    fn compute_part() {
        let mut real_pos = UInt64Chunked::new_vec(
            "real_pos",
            vec![10, 50, 110, 326512443305, 326512443305, 224],
        );
        let mut ref_seq = Utf8Chunked::new("ref", vec!["A", "C", "T", "G", "GA", "AC"]);
        let mut alt_seq = Utf8Chunked::new("alt", vec!["G", "T", "C", "A", "", "CATGAGCGGACTGACC"]);

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
