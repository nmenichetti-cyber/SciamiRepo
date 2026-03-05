#include "FlagClass.h"

bool Flag::IsTriple08(ULong64_t CH){

    return ((CH >> 0) & 1u);

}

bool Flag::IsTriple06(ULong64_t CH){

    return ((CH >> 2) & 1u);

}

bool Flag::IsTriple04(ULong64_t CH){

    return ((CH >> 4) & 1u);

}

bool Flag::IsDouble08(ULong64_t CH){

    return ((CH >> 1) & 1u);

}

bool Flag::IsDouble06(ULong64_t CH){

    return ((CH >> 3) & 1u);

}

bool Flag::IsDouble04(ULong64_t CH){

    return ((CH >> 5) & 1u);

}

bool Flag::IsResetKW(ULong64_t CH){

    if (CH == 2147483648) {return true;}

    else {return false;}

}


