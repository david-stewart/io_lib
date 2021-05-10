#include "ioTags.h"

const char* ioTag_ax(int i) {
    switch (i) {
        case 0: return "#phi_{Trigger}"; break;
        case 1: return "#phi_{Recoil}"; break;
        case 2: return "#phi_{Transverse}"; break;
        case 100: return "#phi_{all}"; break;
        default: return "bad ioTag_ax";
    }
};

const char* ioTag_ea(int i) {
    switch (i) {
        case 0: return "#EA_{lo:70-90%}"; break;
        case 1: return "#EA_{med:30-70%}"; break;
        case 2: return "#EA_{hi:0-30%}"; break;
        case 10: return "#EA_{lo:70-100%}"; break;
        case 20: return "#EA_{lo:90-100%}"; break;
        case 100: return "#EA_{all}"; break;
        default:
                return "bad ioTag_ea";
    }
};

