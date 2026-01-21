#ifndef SVP_H
#define SVP_H

// #include "pipeline.h"
#include <vector>
#include <cstdint>
#include <cassert>
#include <iostream>

#include "debug.h"
 #include "parameters.h"
// #include "parameters.cc"
#include <cmath>
#include "svp_row_t.h"

using namespace std;

// class pipeline_t; // Forward declaration

class svp
{
public:
    // payload *PAY;

    // Member variables
    uint64_t VPQ_size;
    uint64_t index_bits;
    uint64_t tag_bits;
    uint64_t conf_max;
    uint64_t oracleconf;
    int print_count;
    uint64_t perf_val_pred;


    //meas
     uint64_t vpmeas_inelegible;
     uint64_t vpmeas_elegible;
     uint64_t vpmeas_miss;
     uint64_t vpmeas_conf_corr;
     uint64_t vpmeas_conf_incorr;
     uint64_t vpmeas_unconf_corr;
     uint64_t vpmeas_unconf_incorr;
    

    //svp size
    uint64_t svp_size;

    // uint64_t PERF_VAL_PRED;
    // pipeline_t* pipeline;

    // db_t *actual;

    bool stall_v;

    // stats
   // uint64_t vpmeas_inelegible;
   // uint64_t vpmeas_elegible;
    // end stats

    uint64_t head;
    uint64_t tail;
    bool head_phase;
    bool tail_phase;

    struct svp_row_t {
        uint64_t tag;
        uint64_t conf;
        int64_t retired_value;
        int64_t stride;
        uint64_t instance = 0;
        //bool valid = false;
    };


    /*struct svp_row_t {
       uint64_t tag;
       uint64_t conf;
       int64_t retired_value;
       int64_t stride;
       uint64_t instance = 0;
       bool valid = false;
   };*/

    //  // SVP table
    std::vector<svp_row_t> svp_table;

    struct VPQEntry
    {
        int64_t actual_value;
        uint64_t pc;
    };

    // VPQ management
    std::vector<VPQEntry> vpq_table;

    svp(uint64_t vpq_size,
        uint64_t index_bits,
        uint64_t tag_bits,
        uint64_t conf_max,
        uint64_t perf_val_pred,
        uint64_t oracleconf)
    {

        this->VPQ_size = vpq_size;
        this->index_bits = index_bits;
        this->tag_bits = tag_bits;
        // this->tag_bits = 64 - index_bits;
        this->conf_max = conf_max;
        this->oracleconf = oracleconf;
	this->perf_val_pred = perf_val_pred;
        // this->PERF_VAL_PRED = PERF_VAL_PRED;
        // this->pip = pip;
        // this->PAY = PAY;

        stall_v = 0;

        cout << "Constructor called with parameters:" << endl;
        cout << "VPQ_size: " << VPQ_size << endl;
        cout << "index_bits: " << index_bits << endl;
        cout << "tag_bits: " << tag_bits << endl;
        cout << "conf_max: " << conf_max << endl;
        cout << "oracleconf: " << oracleconf << endl;

        int print_count = 0;

        svp_size = (1<<index_bits);

        vpq_table.resize(VPQ_size);
        svp_table.resize(svp_size); // changed to 1 from 2

        head = 0;
        tail = 0;
        head_phase = false;
        tail_phase = false;

        //vpmeas_inelegible = 0;
        //vpmeas_elegible = 0;

	//meas

	vpmeas_inelegible = 0 ;
    	vpmeas_elegible = 0 ;
    	vpmeas_miss = 0;
   	vpmeas_conf_corr = 0;
    	vpmeas_conf_incorr = 0;
    	vpmeas_unconf_corr = 0;
    	vpmeas_unconf_incorr = 0;
        
    }

    int64_t predict(uint64_t PC)
    {

        uint64_t index = mask_n_bits(PC >> 2, index_bits);
        uint64_t tag = mask_n_bits(PC >> (index_bits + 2), tag_bits);
        int64_t prediction;

        assert(svp_table[index].tag == tag);

        svp_table[index].instance = svp_table[index].instance + 1;
        // prediction.valid = true;
        prediction = ((svp_table[index].instance) * (svp_table[index].stride)) + svp_table[index].retired_value;
        // cout << "[PREDICT] TAG: "<< svp_table[index].tag << " instance: " << svp_table[index].instance <<" Prediction: "<<prediction<< endl;
        return prediction;
    }

    uint64_t mask_n_bits(uint64_t value, int n)
    {
        if (n >= 64)
            return value; // All bits kept
        uint64_t mask = (1ULL << n) - 1;
        return value & mask;
    }

    bool stall_vpq(uint64_t bundle_vpq)
    {
        if (vpq_is_full())
        {
            stall_v = 1;
            return true;
        }
        if (head_phase == tail_phase)
        {
            stall_v = ((VPQ_size - (tail - head)) < bundle_vpq);
            return ((VPQ_size - (tail - head)) < bundle_vpq);
        }
        else
        {
            stall_v = ((head - tail) < bundle_vpq);
            return ((head - tail) < bundle_vpq);
        }
    }

    // print_val_meas

   // void print_val_meas(FILE *stats_log)
    //{


	    

        /* fprintf(FP, "\n-VPU MEASUREMENTS-----------------------------------\n");
        fprintf(FP, "vpmeas_ineligible         :\n");

        for (int i = 0; i < 1024; i++)
        {
            fprintf(FP, "%10d: %10lx %6lu %16lu %8ld %10lu\n",
                    i,
                    svp_table[i].tag,
                    svp_table[i].conf,
                    svp_table[i].retired_value,
                    static_cast<int64_t>(svp_table[i].stride), // for negative strides
                    svp_table[i].instance);
        }*/


  //  }



	void print_val_meas(FILE *fp) {
	        uint64_t total = vpmeas_inelegible + vpmeas_elegible;
	    
	        double perc_vpmeas_ineligible = (total == 0) ? 0.0 : (vpmeas_inelegible * 100.0) / total;
	        double perc_vpmeas_eligible   = (total == 0) ? 0.0 : (vpmeas_elegible * 100.0) / total;
	    
	        double perc_miss             = (total == 0) ? 0.0 : (vpmeas_miss * 100.0) / total;
	        double perc_conf_corr        = (total == 0) ? 0.0 : (vpmeas_conf_corr * 100.0) / total;
	        double perc_conf_incorr      = (total == 0) ? 0.0 : (vpmeas_conf_incorr * 100.0) / total;
	        double perc_unconf_corr      = (total == 0) ? 0.0 : (vpmeas_unconf_corr * 100.0) / total;
	        double perc_unconf_incorr    = (total == 0) ? 0.0 : (vpmeas_unconf_incorr * 100.0) / total;
	    
	        fprintf(fp, "VPU MEASUREMENTS-----------------------------------\n");
	        fprintf(fp, "vpmeas_ineligible         : %10lu (%.2f%%) // Not eligible for value prediction.\n", vpmeas_inelegible, perc_vpmeas_ineligible);
	        fprintf(fp, "vpmeas_eligible           : %10lu (%.2f%%) // Eligible for value prediction.\n", vpmeas_elegible, perc_vpmeas_eligible);
	        fprintf(fp, "vpmeas_miss               : %10lu (%.2f%%) // VPU was unable to generate a value prediction (e.g., SVP miss).\n", vpmeas_miss, perc_miss);
	        fprintf(fp, "vpmeas_conf_corr          : %10lu (%.2f%%) // VPU generated a confident and correct value prediction.\n", vpmeas_conf_corr, perc_conf_corr);
	        fprintf(fp, "vpmeas_conf_incorr        : %10lu (%.2f%%) // VPU generated a confident and incorrect value prediction. (MISPREDICTION)\n", vpmeas_conf_incorr, perc_conf_incorr);
	        fprintf(fp, "vpmeas_unconf_corr        : %10lu (%.2f%%) // VPU generated an unconfident and correct value prediction. (LOST OPPORTUNITY)\n", vpmeas_unconf_corr, perc_unconf_corr);
	        fprintf(fp, "vpmeas_unconf_incorr      : %10lu (%.2f%%) // VPU generated an unconfident and incorrect value prediction.\n", vpmeas_unconf_incorr, perc_unconf_incorr);
	    }
   


    void print_svp(FILE *FP, uint64_t cycles)
    {
        fprintf(FP, "\n-------- DEBUG: retired instruction count = %lu\n", cycles);
        fprintf(FP, "SVP entry #:   tag(hex)   conf   retired_value   stride   instance\n");

        for (int i = 0; i < svp_size; i++)
        {
            fprintf(FP, "%10d: %10lx %6lu %16lu %8ld %10lu\n",
                    i,
                    svp_table[i].tag,
                    svp_table[i].conf,
                    svp_table[i].retired_value,
                    static_cast<int64_t>(svp_table[i].stride), // for negative strides
                    svp_table[i].instance);
        }
    }

    void print_vpq(FILE *FP)
    {
        fprintf(FP, "\nVPQ entry #:   PC(hex)   PCtag(hex)   PCindex(hex)\n");

        uint64_t current = head;
        bool current_phase = head_phase;
        int entry_number = 0;

        while (current != tail || current_phase != tail_phase)
        {
            uint64_t pc = vpq_table[current].pc;
            uint64_t pc_index = mask_n_bits(pc >> 2, index_bits);
            uint64_t pc_tag = mask_n_bits(pc >> (index_bits + 2), tag_bits);

            fprintf(FP, "%10d: %10lx %12lx %14lx\n",
                    entry_number,
                    pc,
                    pc_tag,
                    pc_index);

            current = (current + 1) % vpq_table.size();
            if (current == 0)
                current_phase = !current_phase;
            entry_number++;
        }

        // fprintf(FP, "\nInstance counter test: %lu \n\n",count_vpq_entries_with_pc(0x1040c));
    }

    bool vpq_is_empty() const
    {
        return (head == tail) && (head_phase == tail_phase);
    }

    bool vpq_is_full() const
    {
        return (head == tail) && (head_phase != tail_phase);
    }

    VPQEntry vpq_pop()
    {
        if (vpq_is_empty())
        {
            assert(0 && "VPQ pop failed: queue empty");
        }
        VPQEntry entry = vpq_table[head];
        head = (head + 1) % VPQ_size;
        if (head == 0)
            head_phase = !head_phase;
        return entry;
    }

    uint64_t count_vpq_entries_with_pc(uint64_t target_pc)
    {
        uint64_t count = 0;
        uint64_t current = head;
        bool current_phase = head_phase;

        while (current != tail)
        {
            if (vpq_table[current].pc == target_pc)
            {
                count++;
            }
            current = (current + 1) % VPQ_size;
            if (current == 0)
            {
                current_phase = !current_phase;
            }
        }
        return count;
    }

    void train_svp()
    {

        uint64_t current = head;
        bool current_phase = head_phase;

        VPQEntry popped_value = vpq_pop();

        uint64_t PC = popped_value.pc;
        int64_t actual_value = popped_value.actual_value;

        uint64_t svp_index = mask_n_bits(PC >> 2, index_bits);
        uint64_t svp_tag = mask_n_bits(PC >> (index_bits + 2), tag_bits);
        // one head passes a vpq entry we need to train - call this function then
        // to train:
        // 1) update tag
        // 2) check if curr_val - prev_val == stride ---> incr conf    else: conf = 0
        // 3) stride = vpq_val - retired value
        // 4) set retired_value

        if (svp_table[svp_index].tag == svp_tag)
        {
            if ((static_cast<int64_t>(actual_value) - static_cast<int64_t>(svp_table[svp_index].retired_value)) == static_cast<int64_t>(svp_table[svp_index].stride))
            {
                if (svp_table[svp_index].conf < conf_max)
                {
                    svp_table[svp_index].conf++;
                    assert(svp_table[svp_index].conf <= conf_max);
                    // assert(conf_max == 3);
                }
            }
            else
            {   
                svp_table[svp_index].conf = 0; // Figure out if we need to decrement conf or set to 0
                svp_table[svp_index].stride = static_cast<int64_t>(actual_value) - static_cast<int64_t>(svp_table[svp_index].retired_value);
            }

            svp_table[svp_index].retired_value = actual_value;
            // svp_table[svp_index].instance = count_vpq_entries_with_pc(PC);
            if (svp_table[svp_index].instance > 0)
                svp_table[svp_index].instance = count_vpq_entries_with_pc(PC);

        }
        else
        {
            svp_table[svp_index].retired_value = actual_value;
            svp_table[svp_index].tag = svp_tag;
            svp_table[svp_index].conf = 0; // Figure out if we need to decrement conf or set to 0
            svp_table[svp_index].stride = static_cast<int64_t>(actual_value) ;
            svp_table[svp_index].instance = count_vpq_entries_with_pc(PC);
        }
      //  svp_table[svp_index].valid = true;
        // assert(count_vpq_entries_with_pc(PC) == svp_table[svp_index].instance );
    }
    

    bool is_val_pred_in_svp(uint64_t PC)
    {
        // uint64_t PC = pipeline->PAY.buf[index].pc;
        uint64_t svp_index = mask_n_bits(PC >> 2, index_bits);
        uint64_t svp_tag = mask_n_bits(PC >> (index_bits + 2), tag_bits);
        if (svp_table[svp_index].tag == svp_tag )
            return true;
        else
            return false;
    }

    bool is_confident(uint64_t PC, uint64_t actual_value)
    { // actual value for oracleconf
        // uint64_t PC = pipeline->PAY.buf[index].pc;
        uint64_t svp_index = mask_n_bits(PC >> 2, index_bits);
        uint64_t svp_tag = mask_n_bits(PC >> (index_bits + 2), tag_bits);

        if (oracleconf)
        {
            int64_t predicted_value = ((svp_table[svp_index].instance) * (svp_table[svp_index].stride)) + svp_table[svp_index].retired_value;

            if (predicted_value == actual_value)
            {
                return 1;
            }
            else
                return 0;
        }
        else
        {
            //cout<<"in is_confident: svp_table[svp_index].conf: "<<svp_table[svp_index].conf<<" conf_max: "<<conf_max<<endl;
            return (svp_table[svp_index].conf == conf_max);
        }
    }

    //  void create_svp_entry(uint64_t PC , uint64_t instance) {

    //     uint64_t svp_index =  mask_n_bits(PC>>2, index_bits );
    // 	uint64_t svp_tag =  mask_n_bits(PC>>(index_bits+2) , tag_bits );

    //     svp_table[svp_index].tag = svp_tag;
    //     svp_table[svp_index].conf = 0;
    //     svp_table[svp_index].retired_value = 0;
    //     svp_table[svp_index].stride = 0;
    //     svp_table[svp_index].instance  = instance;
    //  }

    uint64_t create_vpq_entry(reg_t PC)
    {
        if (vpq_is_full())
        {
            assert(0 && "VPQ write failed: queue full");
        }
        vpq_table[tail].pc = PC;
        uint64_t vpq_index = tail;
        tail = (tail + 1) % VPQ_size;
        if (tail == 0)
            tail_phase = !tail_phase;
        return vpq_index;
    }

    void update_vpq_entry(uint64_t vpq_index, int64_t ret_val)
    { // to add executed value to the corresponding entry

        // write assertion to check if this index is btw head and tail
        // assert(vpq_index >= head && "vpq entry to update is not a valid entry" );
        // assert(vpq_index < tail && "vpq entry to update is not a valid entry" );

        vpq_table[vpq_index].actual_value = ret_val;
    }

    // recovery

    void vpq_partial_rollback(uint64_t vpq_tail_cp, bool vpq_tail_phase_cp)
    {
        if((vpq_tail_cp == head) && (vpq_tail_phase_cp == head_phase)){
            vpq_full_rollback();
        }
        else{
            uint64_t local_tail = vpq_tail_cp;
            uint64_t index;
            bool in_svp;
            // we forward walk the VPQ
            // checkpointed tail is bullshit entry
            // any entry following it also bullshit
            // our tail doesnt point to bullshit it just points to empty so do not use it to remove bullshit from svp
            while (local_tail != tail)
            { // when == tail we can stop undoing
                index = mask_n_bits(vpq_table[local_tail].pc >> 2, index_bits);
                in_svp = is_val_pred_in_svp(vpq_table[local_tail].pc);

                if (in_svp)
                { // if count in svp exists - it incremented when it when predicted it so decrement it
                    if (svp_table[index].instance > 0){
                        svp_table[index].instance--;
                    }
                }
                local_tail = (local_tail + 1) % VPQ_size;
            }
            tail_phase = vpq_tail_phase_cp;
            tail = vpq_tail_cp; // were back to where we were
        }
    }

    void vpq_full_rollback()
    {
        tail = head;
        tail_phase = head_phase;

        for (int i = 0; i < svp_size; i++)
        {
            svp_table[i].instance = 0;
        }
    }
};

#endif // SVP_H
