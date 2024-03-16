#ifndef _SIMPLEX_INL_
#define _SIMPLEX_INL_

inline bool Simplex::is_face(const Simplex &other) const
{
    if ((bloom_filter_ & other.bloom_filter()) != bloom_filter_)
    {
        return false;
    }
    return std::includes(other.vertices().begin(), other.vertices().end(), vertices_.begin(), vertices_.end());
}

inline std::bitset<BLOOM_FILTER_SIZE> Simplex::bloom_filter() const
{
    return bloom_filter_;
}

#endif