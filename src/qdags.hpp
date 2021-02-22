#ifndef INCLUDED_QDAGS
#define INCLUDED_QDAGS

#include <sdsl/bit_vectors.hpp>
#include "se_quadtree.hpp"

#include<bits/stdc++.h>
#include<ratio>
#include<chrono>
#include<ctime>

    extern high_resolution_clock::time_point start_rank, stop_rank;
    extern double total_time_rank;       
    extern duration<double> time_span_rank;

typedef uint8_t type_mapping_M;

bool compare_pairs(const pair<uint64_t, uint64_t>&i, const pair<uint64_t, uint64_t>&j)
{
    return i.second < j.second;
}


class qdag 
{
    public:
        
        typedef vector<uint64_t> att_set;

	 private:
        se_quadtree* Q;
	            
        type_mapping_M*    M;    // mapping
        
        att_set      attribute_set;
        
        uint64_t     grid_side;
        
        uint16_t     Msize;  // number of children of every qdag node

        bool         is_extended_qdag;
        
        vector<vector<type_mapping_M>*> M_prime;

        int32_t tab_extend_5[16];   // queries of 5 attributes, i.e., dimension 2^5=32
        int32_t tab_extend_4[16];
        int32_t tab_extend_3[16]; 

    public:
        
        qdag() = default;

        uint64_t size()
	{
	    uint64_t s = Q->size() + Msize*sizeof(uint16_t) + attribute_set.size()*sizeof(uint64_t) 
                        + M_prime.size()*sizeof(vector<type_mapping_M>*) + sizeof(uint64_t);

            //for (uint64_t i = 0; i < M_prime.size(); i++)
            //    s += M_prime[i]->size()*sizeof(type_mapping_M);

	    return s;
        }
      
       
        void setAtts(uint64_t att1, uint64_t att2)
	{
	    attribute_set[0] = att1;
	    attribute_set[1] = att2;
	}

        qdag (const qdag &_Q)
	{
	    this->Q = _Q.Q;
            _Q.Q->inc_ref_count();
            this->M = _Q.M;
	    for (uint64_t i=0; i < _Q.attribute_set.size(); i++)
	        this->attribute_set.push_back(_Q.attribute_set[i]);

            this->grid_side = _Q.grid_side;
	    this->Msize = _Q.Msize;
	    this->is_extended_qdag = _Q.is_extended_qdag;

	}

        qdag(std::vector<std::vector<uint64_t>> &points, 
             att_set &_attribute_set, 
             const uint64_t _grid_side,
             uint8_t k, uint8_t d
            )
        {            

           Msize = std::pow(k, d);           
           
           M = new type_mapping_M[Msize];
           
           uint64_t i, j;
           
           for (i = 0; i < Msize; i++)
               M[i] = i;  // identity mapping 
           
           attribute_set = _attribute_set;
           
           vector<uint64_t> tuple_aux(d);
           vector<pair<uint64_t,uint64_t>> map_sort_att(d);
           
           for (i = 0; i < d; i++)
               map_sort_att[i] = make_pair(i, attribute_set[i]);           
           
           std::sort(map_sort_att.begin(), map_sort_att.end(), compare_pairs);
           
           for (i = 0; i < points.size(); i++) {
               //if (i%1000000==0) cout << i << endl;  
               for (j = 0; j < d; j++)
                   tuple_aux[j] = points[i][map_sort_att[j].first]; 
               
               for (j = 0; j < d; j++)
                   points[i][j] = tuple_aux[j];
                  
           }

           std::sort(attribute_set.begin(), attribute_set.end());            

           //cout << "Construyendo el quadtree" << endl;
           Q = new se_quadtree(points, _grid_side, k, d);          
           
           grid_side = _grid_side;
           is_extended_qdag = false;
           
           //M_prime.reserve(Msize);

           //for (uint64_t i = 0; i < Msize; i++)              
           //    M_prime.push_back(new std::vector<type_mapping_M>());
           
           //for (uint64_t i = 0; i < Msize; i++)
           //    M_prime[i]->push_back(i);            
           
        }
      

        qdag(vector<uint64_t> bv[], 
             att_set &_attribute_set, 
             const uint64_t _grid_side,
             uint8_t k, uint8_t d
				)
        {
           Q = new se_quadtree(bv, _grid_side, k, d);

           Msize = std::pow(k, d);           
           
           M = new type_mapping_M[Msize];
           
           for (uint64_t i = 0; i < Msize; i++)
               M[i] = i;  // identity mapping 
           
           attribute_set = _attribute_set;
           //std::sort(attribute_set.begin(), attribute_set.end());            
           grid_side = _grid_side;
           is_extended_qdag = false;
 
           //M_prime.reserve(Msize);
           
           //for (uint64_t i = 0; i < Msize; i++)              
           //    M_prime.push_back(new std::vector<type_mapping_M>());
           
           //for (uint64_t i = 0; i < Msize; i++) {               
           //    M_prime[M[i]]->push_back(i);
           //}
                   	
        }            


        /*qdag(qdag &q, att_set &_attribute_set)
        {
            this->Q = q.Q;
	    Msize = q.Msize;
	    M = new type_mapping_M[Msize];
	    for (uint64_t i = 0; i < Msize; i++)
	        M[i] = q.M[i];

            attribute_set = _attribute_set;
	    std::sort(attribute_set.begin(), attribute_set.end());
	}*/
      
        ~qdag() 
        {
             //if (Q && !is_extended_qdag) {
             //    delete Q;  
             //    Q = NULL;
             //}             
             if (is_extended_qdag) delete M;       
        }
       

        qdag* extend(att_set &attribute_set_A)
        {
            uint16_t dim = attribute_set_A.size();
            uint16_t dim_prime = attribute_set.size();
            uint64_t p = std::pow(Q->getK(), dim);            
            
            type_mapping_M* _M = new type_mapping_M[p];
            
            uint64_t mask;
            uint64_t i, i_prime;            

            for (i = 0; i < p; ++i) {
                mask = 1<<(dim_prime-1);
                i_prime = 0;
               
                for (uint16_t j = 0; j < dim_prime; ++j) {
                    if (i & (1 << (dim-attribute_set[j]-1)))
                        i_prime |= mask;
                
                    mask >>= 1;
                }
            
                _M[i] = M[i_prime];
            }
             
            qdag* q = new qdag();
            
            q->Q = this->Q;
            q->M = _M;
            q->attribute_set = attribute_set_A;
            q->grid_side = this->grid_side; 
            q->is_extended_qdag = true;
            q->Msize = p; // this.Msize;
            
            return q;
        }


        uint64_t nAttr() 
        {
            return attribute_set.size();        
        }


        uint64_t getAttr(uint64_t i) 
        {
            return attribute_set[i];        
        }
        
        
        uint64_t getGridSide() 
        {
            return grid_side;        
        }

        
        uint64_t getHeight() 
        {
            return Q->getHeight();        
        }


        uint8_t getK() 
        {
            return Q->getK();     
        }
        
        
        uint16_t nChildren()
        {
            return Msize;        
        }
     
              
        uint64_t getKD()
        {
            return Q->getKD();
        }


        uint16_t getM(uint16_t i)
        {
            return M[i];        
        }


        // This is for a binary relation, i.e., a k^2-tree with 4 children per node
        void createTableExtend5()
        {
            uint64_t i, j;
            uint32_t x, B;

            if (Q->getKD() == 2)
                B = 4;
            else
                B = 16;
            
            for (i = 0; i < B; i++) {
                x = 0;
                for(j = 0; j < 32; j++)            
                    if (i& (1<<M[j]))
                       x = (x << 1) | 1;
                    else
                       x = (x << 1);
                
                tab_extend_5[i] = x;
                //cout << x << endl;
            }                
        }
        
        // This is for a binary relation, i.e., a k^2-tree with 4 children per node
        void createTableExtend4()
        {
            uint64_t i, j;
            uint32_t x, B;

            if (Q->getKD() == 2)
                B = 4;
            else
                B = 16;

            for (i = 0; i < B; i++) {
                x = 0;
                for(j = 0; j < 16; j++)
                    if (i& (1<<M[j])) 
                       x = (x << 1) | 1;
                    else 
                       x = (x << 1);

                tab_extend_4[i] = x<<16;
                //cout << std::hex << x << endl;
            }
        }

        // This is for a binary relation, i.e., a k^2-tree with 4 children per node
        void createTableExtend3()
        {
            uint64_t i, j;
            uint32_t x, B=16;

            //if (Q->getKD() == 2) 
            //    B = 4;
            //else
            //    B = 16;
 
            for (i = 0; i < B; i++) {
                x = 0;
                for(j = 0; j < 8; j++)
                    if (i& (1<<M[j]))
                       x = (x << 1) | 1;
                    else
                       x = (x << 1);

                tab_extend_3[i] = x<<24;
                //cout << std::hex << x << endl;
            }
        }

        inline uint32_t materialize_node_3(uint64_t level, uint64_t node, uint64_t* rank_vector) {
            uint64_t r = Q->rank(level, node);
            return tab_extend_3[Q->get_node(level, node, rank_vector, r)];
        }


        inline uint32_t materialize_node_4(uint64_t level, uint64_t node, uint64_t* rank_vector) {
            uint64_t r = Q->rank(level, node);
            return tab_extend_4[Q->get_node(level, node, rank_vector, r)];
        }


        inline uint32_t materialize_node_5(uint64_t level, uint64_t node, uint64_t* rank_vector) {
            uint64_t r = Q->rank(level, node);            
            return tab_extend_5[Q->get_node(level, node, rank_vector, r)]; 
        }


        inline uint32_t materialize_node_3_lastlevel(uint64_t level, uint64_t node) {
            return tab_extend_3[Q->get_node_lastlevel(level, node)];
        }


        inline uint32_t materialize_node_4_lastlevel(uint64_t level, uint64_t node) {
            return tab_extend_4[Q->get_node_lastlevel(level, node)];
        }
       

        inline uint32_t materialize_node_5_lastlevel(uint64_t level, uint64_t node) {
            return tab_extend_5[Q->get_node_lastlevel(level, node)];
        }
 
           
        //void print(std::ofstream &ofs)
        //{
        //    Q->print(ofs);        
        //}
};

#endif
