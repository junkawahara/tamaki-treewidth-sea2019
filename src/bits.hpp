// Dynamic-size bitset used for vertex sets.
#pragma once
#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

class Bits {
 public:
  std::vector<uint64_t> w;

  Bits() = default;
  // n >= 0.  The word count is computed in size_t so that n near INT_MAX does not overflow.
  explicit Bits(int n) : w((static_cast<size_t>(n) + 63) >> 6, 0) {}

  void set(int i) { w[i >> 6] |= 1ULL << (i & 63); }
  void reset(int i) { w[i >> 6] &= ~(1ULL << (i & 63)); }
  bool test(int i) const { return (w[i >> 6] >> (i & 63)) & 1ULL; }

  int count() const {
    int c = 0;
    for (uint64_t x : w) c += std::popcount(x);
    return c;
  }
  bool none() const {
    for (uint64_t x : w)
      if (x) return false;
    return true;
  }
  bool any() const { return !none(); }

  // Smallest element, or -1 when empty.
  int first() const {
    for (size_t i = 0; i < w.size(); ++i)
      if (w[i]) return static_cast<int>(i * 64 + std::countr_zero(w[i]));
    return -1;
  }
  // Smallest element of (*this & o), or -1.
  int firstCommon(const Bits& o) const {
    for (size_t i = 0; i < w.size(); ++i) {
      uint64_t x = w[i] & o.w[i];
      if (x) return static_cast<int>(i * 64 + std::countr_zero(x));
    }
    return -1;
  }

  template <class F>
  void forEach(F f) const {
    for (size_t i = 0; i < w.size(); ++i) {
      uint64_t x = w[i];
      while (x) {
        f(static_cast<int>(i * 64 + std::countr_zero(x)));
        x &= x - 1;
      }
    }
  }
  std::vector<int> toVector() const {
    std::vector<int> v;
    forEach([&](int i) { v.push_back(i); });
    return v;
  }

  Bits& operator|=(const Bits& o) {
    for (size_t i = 0; i < w.size(); ++i) w[i] |= o.w[i];
    return *this;
  }
  Bits& operator&=(const Bits& o) {
    for (size_t i = 0; i < w.size(); ++i) w[i] &= o.w[i];
    return *this;
  }
  // Set difference.
  Bits& operator-=(const Bits& o) {
    for (size_t i = 0; i < w.size(); ++i) w[i] &= ~o.w[i];
    return *this;
  }

  bool subsetOf(const Bits& o) const {
    for (size_t i = 0; i < w.size(); ++i)
      if (w[i] & ~o.w[i]) return false;
    return true;
  }
  bool intersects(const Bits& o) const {
    for (size_t i = 0; i < w.size(); ++i)
      if (w[i] & o.w[i]) return true;
    return false;
  }
  bool operator==(const Bits& o) const { return w == o.w; }
  bool operator!=(const Bits& o) const { return w != o.w; }

  size_t hash() const {
    uint64_t h = 0x9E3779B97F4A7C15ULL;
    for (uint64_t x : w) {
      h ^= x + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
      h *= 0xBF58476D1CE4E5B9ULL;
      h ^= h >> 31;
    }
    return static_cast<size_t>(h);
  }
};

inline Bits operator|(Bits a, const Bits& b) { return a |= b; }
inline Bits operator&(Bits a, const Bits& b) { return a &= b; }
inline Bits operator-(Bits a, const Bits& b) { return a -= b; }

inline int unionCount(const Bits& a, const Bits& b) {
  int c = 0;
  for (size_t i = 0; i < a.w.size(); ++i) c += std::popcount(a.w[i] | b.w[i]);
  return c;
}

struct BitsHash {
  size_t operator()(const Bits& b) const { return b.hash(); }
};
