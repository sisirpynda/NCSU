#ifndef SVP_ROW_H
#define SVP_ROW_H

struct svp_row_t {
        uint64_t tag;
        uint64_t conf;
        int64_t retired_value;
        int64_t stride;
        uint64_t instance = 0;
        //bool valid = false;
    };

#endif // SVP_ROW_H
