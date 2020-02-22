
#include "Ring/BPASRing.hpp"

RingProperties::RingProperties() {
    prop = 0;
}

RingProperties::RingProperties(RingProperty p) {
    // n = 1; 
    // prop = new RingProperty[n];
    // props[0] = p;
    prop = p;
}

RingProperties::RingProperties(std::vector<RingProperty> v) {
    prop = 0;
    for (int i = 0 ; i < v.size(); ++i) {
        prop |= v[i];
    }
}

inline bool RingProperties::has(RingProperty p) {
    return prop & p; 
    // for (int i = 0; i < n; ++i) {
    //     if (prop[i] & p) {
    //         return 1;
    //     }
    // }
    // return 0;
}

inline bool RingProperties::has(const RingProperties& p) {
    // for (int i = 0; i < p.n; ++i) {
    //     if (!this.has(p.prop[i])) {
    //         return 0;
    //     }
    // }
    // return 1;
    return prop & p.prop;
}
