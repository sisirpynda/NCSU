// svp.cc
#include "svp.h"
//#include "pipeline.h"
#include <cmath>




// uint64_t svp::predict(db_t &actual , payload *PAY, unsigned int index ) {
//     if (PERF_VAL_PRED ) {
//         //db_t* actual = pip->get_pipe()->peek(PAY->buf[index].db_index);
//         uint64_t predicted_val = actual.a_rdst[0].value;
        
//         // Write to VPQ
//         if (!vpq_is_full()) {
//             VPQEntry entry;
//             entry.actual_value = predicted_val;
//             entry.pc = PAY->buf[index].pc;
//             vpq_table[tail] = entry;
//             tail = (tail + 1) % VPQsize;
//             if (tail == 0) tail_phase = !tail_phase;
//         }
//         return predicted_val;
//     }
//     return 0;
// }

// void svp::train_svp(unsigned int index, uint64_t vpq_val) {
//     uint64_t PC = PAY.buf[index].pc;
//     uint64_t svp_index = mask_n_bits(PC, index_bits);
//     uint64_t svp_tag = mask_n_bits(PC >> index_bits, tag_bits);

//     svp_table[svp_index].tag = svp_tag;
//     if (((vpq_val - svp_table[svp_index].retired_value) == svp_table[svp_index].stride) && 
//         (svp_table[svp_index].conf != confmax)) {
//         svp_table[svp_index].conf++;
//     } else if (svp_table[svp_index].conf > 0) {
//         svp_table[svp_index].conf--;
//     }
//     svp_table[svp_index].stride = vpq_val - svp_table[svp_index].retired_value;
//     svp_table[svp_index].retired_value = vpq_val;
// }

// bool svp::is_val_pred_in_svp(unsigned int index) {
//     uint64_t PC = PAY.buf[index].pc;
//     uint64_t svp_index = mask_n_bits(PC, index_bits);
//     uint64_t svp_tag = mask_n_bits(PC >> index_bits, tag_bits);
//     return svp_table[svp_index].tag == svp_tag;
// }

// uint64_t svp::mask_n_bits(uint64_t value, int n) {
//     if (n >= 64) return value;
//     uint64_t mask = (1ULL << n) - 1;
//     return value & mask;
// }

// // Remaining VPQ management functions
// void svp::create_vpq_entry(uint64_t PC) {
//     if (!vpq_is_full()) {
//         VPQEntry entry;
//         entry.pc = PC;
//         vpq_table[tail] = entry;
//         tail = (tail + 1) % VPQsize;
//         if (tail == 0) tail_phase = !tail_phase;
//     }
// }

// svp::VPQEntry svp::vpq_pop() {
//     assert(!vpq_is_empty());
//     VPQEntry entry = vpq_table[(int)head];
//     head = (head + 1) % VPQsize;
//     if (head == 0) head_phase = !head_phase;
//     return entry;
// }

// bool svp::vpq_is_empty() const {
//     return (head == tail) && (head_phase == tail_phase);
// }

// bool svp::vpq_is_full() const {
//     return (head == tail) && (head_phase != tail_phase);
// }
