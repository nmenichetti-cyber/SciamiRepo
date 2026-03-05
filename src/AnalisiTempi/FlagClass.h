#ifndef FLAG_H
#define FLAG_H

#include "Rtypes.h"   // per ULong64_t (ROOT)

class Flag {
public:
    bool IsTriple08(ULong64_t CH);
    bool IsTriple06(ULong64_t CH);
    bool IsTriple04(ULong64_t CH);
    bool IsDouble08(ULong64_t CH);
    bool IsDouble06(ULong64_t CH);
    bool IsDouble04(ULong64_t CH);
    bool IsResetKW(ULong64_t CH);
};

#endif