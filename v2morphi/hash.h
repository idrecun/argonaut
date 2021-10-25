#ifndef HASH_H
#define HASH_H

#include <cstdint>

namespace morphi::hash {

uint32_t hash32(uint32_t value) {
    value++;
    value ^= value >> 17;
    value *= 0xed5ad4bb;
    value ^= value >> 11;
    value *= 0xac4c1b51;
    value ^= value >> 15;
    value *= 0x31848bab;
    value ^= value >> 14;
    return value;
}

uint32_t sequential32(uint32_t prev_hash, uint32_t value) {
    return hash32(prev_hash ^ value);
}

void sequential32u(uint32_t& prev_hash, uint32_t value) {
    prev_hash = sequential32(prev_hash, value);
}

uint32_t multiset32(uint32_t prev_hash, uint32_t value) {
    return prev_hash + hash32(value) + 1;
}

void multiset32u(uint32_t& prev_hash, uint32_t value) {
    prev_hash = multiset32(prev_hash, value);
}

} // namespace

#endif // HASH_H
