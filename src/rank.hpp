#ifndef RANK_BV
#define RANK_BV

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector_buffer.hpp>

using namespace sdsl;
using namespace std;


class rank_bv_64
{
    uint64_t* seq;
    uint32_t* block;
    uint64_t u;  //bit vector length
    uint64_t n; // # ones
    
   public:
    rank_bv_64() = default;
    
    rank_bv_64(bit_vector &bv)
    {
        uint64_t i;
        uint8_t byte_mask;
        uint32_t cur_word = 0, count = 0;

        u = bv.size();   
             
        seq = new uint64_t[(u+63)/64]();     
        block = new uint32_t[(u+63)/64]();
        
        for (i = 0; i < u; ++i) {
 
            if (i%64 == 0)
                block[cur_word++] = count; 
                
            if (bv[i]) {
                count++;
                seq[i/64] |= (1L<<(i%64));
            }
            else 
                seq[i/64] &= ~(1L<<(i%64));
            	
        }
        n = count;
    }


    inline uint64_t rank(uint64_t i) 
    {
        return block[i>>6] + bits::cnt(seq[i>>6] & ~(0xffffffffffffffff << (i&0x3f)));  
    }

    inline uint8_t get_4_bits(uint64_t start_pos)
    {
        return ((seq[start_pos >> 6] >>(start_pos & 0x3f) ) & 0x0f);
    }
   
    inline uint8_t get_2_bits(uint64_t start_pos)
    {
        return ((seq[start_pos >> 6] >>(start_pos & 0x3f) ) & 0x03);
    }
 
    // number of bits in the bv
    inline uint64_t size()
    {
        return u;    
    }

    inline uint64_t n_ones()
    {
        return n;
    }

    inline uint64_t size_in_bytes()
    {
        return sizeof(uint64_t)*((u+63)/64) + sizeof(uint32_t)*(u+63)/64 
	       + sizeof(uint64_t*) + sizeof(uint32_t*)
	       + 2*sizeof(uint64_t);
    }

};

#endif
