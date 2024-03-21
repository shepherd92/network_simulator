#ifndef _SIMPLEX_INL_
#define _SIMPLEX_INL_

inline bool Simplex::is_face(const Simplex &other) const
{
    const auto hash_value{hash()};
    if ((hash_value & other.hash()) != hash_value)
    {
        return false;
    }
    return std::includes(other.vertices().begin(), other.vertices().end(), vertices_.begin(), vertices_.end());
}

inline uint64_t Simplex::hash() const
{
    return hash_;
}

#endif