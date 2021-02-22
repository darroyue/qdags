/* 
    Copyright (C) 2020 Diego Arroyuelo
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file se_quadtree.hpp
    \brief se_quadtree.hpp contains a space-efficient quadtree implementation.
    \author Diego Arroyuelo
    \Based on Francisco Montoto's k2-tree implementation (sdsl)

*/
#ifndef SE_QUADTREE
#define SE_QUADTREE

#include <deque>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <fstream> 
#include <sdsl/bit_vectors.hpp>
#include <sdsl/k2_tree_helper.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include "rank.hpp"

//! A k^2-tree
/*! A k^2-tree is a compact tree structure to represent a web graph. The
 *  structure takes advantage of large empty areas of the adjacency matrix of
 *  the graph.
 *
 *  \par References
 *      [1] Brisaboa, N. R., Ladra, S., & Navarro, G. (2009, August):
 *          k2-trees for compact web graph representation. In International
 *          Symposium on String Processing and Information Retrieval
 *          (pp. 18-30). Springer Berlin Heidelberg.
 */

using namespace sdsl;
using namespace std;

class se_quadtree
{
    public:
        typedef k2_tree_ns::idx_type idx_type;
        typedef k2_tree_ns::size_type size_type;

    private:
        rank_bv_64   *bv;  
        
        //rank_support_v<> *bv_rank;

        uint16_t    height;   // number of levels of the tree

        uint8_t     k;
        uint8_t     d;
        size_type   k_d;
        vector<uint64_t> total_ones;
        uint64_t    ref_count;

    protected:
    
 /*! Get the chunk index ([0, k^d[) of a submatrix point.
 *
 * Gets a point in the global matrix and returns its corresponding chunk
 * in the submatrix specified.
 *
 * \param point vector representing the d dimensional point.
 * \param offset vector with the upper-left point of the global matrix
 * \param l size of the chunk at the submatrix.
 * \param k the k parameter from the k^2 tree.
 * \returns the index of the chunk containing the point at the submatrix.
 */    
        uint16_t get_chunk_idx(vector<idx_type> &point, idx_type* offset,
                               size_type l, uint8_t _k)
        {
            uint16_t i, _d = point.size();
            uint16_t total_sum = 0, _k_aux = 1;
                        
            for (i = _d-1; i > 0; --i) {
                total_sum += ((point[i]-offset[i])/l)*_k_aux;
                _k_aux *= _k;
            }
            
            total_sum += ((point[0]-offset[0])/l)*_k_aux;
            return total_sum;
        }
        
 //! Build a space efficient quadtree from a set of d-dimensional points
        /*! This method takes a vector of d-dimensional points 
         *  and the grid-side size. It takes linear time over the amount of
         *  points to build the quadtree representation.
         *  \param edges A non-empty vector with the d-dimensional points
         *  \param size side of the grid containing the points, all the point 
                        coordinates in vector edges must be within [0, size[
         */
        void build_from_edges(std::vector<std::vector<idx_type>>& edges,
                              const size_type size, uint8_t __k, uint8_t __d)
	{

            typedef std::tuple<idx_type, idx_type, size_type, idx_type*> t_part_tuple;
           
            k = __k;
            d = __d;  
            height = std::ceil(std::log(size)/std::log(k));
            height = height > 1 ? height : 1; // If size == 0
            
            bv = new rank_bv_64[height];
            total_ones.reserve(height);
                        
            k_d = std::pow(k, d);
                        
            bit_vector k_t_ = bit_vector(k_d, 0);  // OJO, cuidado con esto

            std::queue<t_part_tuple> q;
            idx_type t = 0, last_level = 0;
            idx_type i, j, r_0, c_0, it, c, r, z;
            size_type l = std::pow(k, height - 1);
            
            std::vector<idx_type> pos_by_chunk(k_d + 1, 0);
            
            idx_type* top_left_point = new idx_type[d]();
            idx_type* top_left_point_aux;            
            
            q.push(t_part_tuple(0, edges.size(), l, top_left_point));

            size_type cur_l = l, cur_level = 0, n_ones = 0;
            
            while (!q.empty()) {
                
                std::vector<idx_type> amount_by_chunk(k_d, 0);
                std::tie(i, j, l, top_left_point) = q.front();
                q.pop();
                
                if (l != cur_l) {                    
                    cur_l = l;
                    k_t_.resize(t);
                    bv[cur_level] = rank_bv_64(k_t_);
                    total_ones[cur_level] =  bv[cur_level].n_ones();    
                    cur_level++;
                    t = 0;
                    k_t_.resize(0);
                    k_t_ = bit_vector(k_d*n_ones, 0);                    
                    n_ones = 0;
                }
                
                // Get size for each chunk
                for (it = i; it < j; it++)              
                    amount_by_chunk[get_chunk_idx(edges[it], top_left_point, l, k)] += 1;
              
                if (l == 1) {
                    for (it = 0; it < k_d; it++,t++)
                        if (amount_by_chunk[it] != 0)
                             k_t_[t] = 1;
                         else
                             k_t_[t] = 0;
                    // At l == 1 no new elements are enqueued 
                    continue;
                }

                // Set starting position in the vector for each chunk
                pos_by_chunk[0] = i;
                
                for (it = 1; it < k_d; it++)
                    pos_by_chunk[it] = pos_by_chunk[it - 1] + amount_by_chunk[it - 1];
                // To handle the last case when it = k_d - 1
                pos_by_chunk[k_d] = j;
                
                // Push to the queue every non zero elements chunk
                for (it = 0; it < k_d; it++,t++)
                    // If not empty chunk, set bit to 1
                    if (amount_by_chunk[it] != 0) {
                        uint16_t p = std::pow(k, d-1);
                        idx_type it_aux = it;                        
                        c = it % k;
                        k_t_[t] = 1;
                        n_ones++;
                        
                        top_left_point_aux = new idx_type[d]();                        
                                               
                        for (z = 0; z < d-1; ++z) {
                            r = it_aux / p;
                            it_aux -= r*p;
                            p /= k;
                            top_left_point_aux[z] = top_left_point[z] + r*l;
                        }
                        
                        top_left_point_aux[d-1] = top_left_point[d-1] + c*l;

                        q.push(t_part_tuple(pos_by_chunk[it],
                                            pos_by_chunk[it + 1],
                                            l/k,
                                            top_left_point_aux));
                    }
                    else
                        k_t_[t] = 0;
                    
                idx_type chunk;
                
                // Sort edges' vector
                for (unsigned ch = 0; ch < k_d; ch++) {
                    idx_type be = ch == 0 ? i : pos_by_chunk[ch - 1];
                    for (it = pos_by_chunk[ch]; it < be + amount_by_chunk[ch];) {
                        chunk = get_chunk_idx(edges[it],top_left_point, l, k);

                        if (pos_by_chunk[chunk] != it)
                            std::iter_swap(edges.begin() + it,
                                           edges.begin() + pos_by_chunk[chunk]);
                        else
                            it++;
                        pos_by_chunk[chunk]++;
                    }
                }
                delete [] top_left_point; 

            }

            k_t_.resize(t);
            bv[height-1] = rank_bv_64(k_t_);
            
            total_ones[height-1] = bv[height-1].n_ones();    
        
        }   

   public:

        se_quadtree() = default;
        
        uint64_t size()
        {
	    uint64_t i, s = 0;
	    for (i = 0; i < height; i++)
	        s += bv[i].size_in_bytes(); 
	    
	    return s + total_ones.size()*sizeof(uint64_t); 
	}


        //! Constructor
        /*! This constructor takes a vector of edges describing the graph
         *  and the graph size. It takes linear time over the amount of
         *  edges to build the k_d representation.
         *  \param edges A vector with all the edges of the graph, it can
         *               not be empty.
         *  \param size grid side, all the point coordinates in edges
                        vector must be within 0 and size ([0, size[).
         */ 
 
        se_quadtree(std::vector<std::vector<idx_type>>& edges,
                const size_type grid_side, uint8_t __k, uint8_t __d)
        {
            assert(grid_side > 0);
            assert(edges.size() > 0);

            build_from_edges(edges, grid_side, __k, __d);
            ref_count = 1;
        }
        
        
        se_quadtree(vector<uint64_t> _bv[], const size_type grid_side, uint8_t _k, uint8_t _d)
        {
            k = _k;
            d = _d;  
            height = std::ceil(std::log(grid_side)/std::log(k));
            height = height > 1 ? height : 1; // If size == 0
            
            bv = new rank_bv_64[height];
            total_ones.reserve(height);
                        
            k_d = std::pow(k, d);
            
            // cout << /*"Join has " <<*/ _bv[height-1].size() << " \t"/* << endl*/;            
            
	    // OJO con lo que sigue, tengo que hacer este constructor de buena manera
		
            /*for (uint64_t i = 0; i < height; i++) {
                bv[i] = new sdsl::sd_vector<>(_bv[i].begin(),_bv[i].end());
                bv_rank[i] = sd_vector<>::rank_1_type(bv[i]);
                bv_select[i] = sd_vector<>::select_1_type(bv[i]);
                total_ones[i] = bv_rank[i].rank(bv[i]->size());    
            }*/
            
        }
        
        
        ~se_quadtree()
        {
            ref_count--;
            if (ref_count == 0)
                delete [] bv;
        }
        
        void inc_ref_count()
        {
	    ref_count++;
	}
 
	
        uint8_t getK()         
        {
            return k;
        }         
        

        uint8_t getD()         
        {
            return d;
        }
        
        
        inline uint64_t rank(uint16_t level, uint64_t node)
        {
            return bv[level].rank(node);
        }

        inline uint8_t get_node_lastlevel(uint16_t level, uint64_t node)
        {
            if (k_d == 4)
                return bv[level].get_4_bits(node);
            else
                return bv[level].get_2_bits(node);
        }

        inline uint8_t get_node(uint16_t level, uint64_t node, uint64_t* rank_array, uint64_t r) 
        {
            uint8_t nd;// = bv[level].get_4_bits(node);
            if (k_d == 4) {
            nd = bv[level].get_4_bits(node);
            switch(nd) {
	        case 0:   break;
		case 1:   rank_array[0] = r+1;
			  break;
		case 2:   rank_array[1] = r+1;
			  break;
		case 3:   rank_array[0] = r+1;
			  rank_array[1] = r+2;
			  break;
		case 4:   rank_array[2] = r+1;
			  break;
		case 5:   rank_array[0] = r+1;
			  rank_array[2] = r+2;
			  break;
		case 6:   rank_array[1] = r+1;
			  rank_array[2] = r+2;
			  break;
		case 7:   rank_array[0] = r+1;
			  rank_array[1] = r+2;
			  rank_array[2] = r+3;
			  break;
		case 8:   rank_array[3] = r+1;
			  break;
		case 9:   rank_array[0] = r+1;
                          rank_array[3] = r+2;
			  break;				    
		case 10:  rank_array[1] = r+1;
			  rank_array[3] = r+2;
			  break;
                case 11:  rank_array[0] = r+1;
			  rank_array[1] = r+2;
			  rank_array[3] = r+3;
			  break;
		case 12:  rank_array[2] = r+1;
			  rank_array[3] = r+2;
			  break;
		case 13:  rank_array[0] = r+1;
			  rank_array[2] = r+2;
			  rank_array[3] = r+3;
			  break;
		case 14:  rank_array[1] = r+1;
			  rank_array[2] = r+2;
			  rank_array[3] = r+3;
			  break;
		case 15:  rank_array[0] = r+1;
			  rank_array[1] = r+2;
			  rank_array[2] = r+3;
			  rank_array[3] = r+4;
                          break;
	    }
            }
            else {
                nd = bv[level].get_2_bits(node);
                switch(nd) {
                    case 0: break;
                    case 1: rank_array[0] = r+1;
                            break;
                    case 2: rank_array[1] = r+1;
                            break;
                    case 3: rank_array[0] = r+1;
                            rank_array[1] = r+2;
                            break;
                }
            }

            return nd;
        }


        void get_children(uint16_t level, uint64_t node, uint64_t* children_array, uint64_t &n_children)
        {
            if (k_d == 4)
            switch(bv[level].get_4_bits(node)) {
	        case 0:   n_children = 0;
		          break;
		case 1:   n_children = 1;
		          children_array[0] = 0;
			  break;
		case 2:   n_children = 1;
		          children_array[0] = 1;
			  break;
		case 3:   n_children = 2;
		          children_array[0] = 0;
			  children_array[1] = 1;
			  break;
		case 4:   n_children = 1;
		          children_array[0] = 2;
			  break;
		case 5:   n_children = 2;
		          children_array[0] = 0;
			  children_array[1] = 2;
			  break;
		case 6:   n_children = 2;
		          children_array[0] = 1;
			  children_array[1] = 2;
			  break;
		case 7:   n_children = 3;
		          children_array[0] = 0;
			  children_array[1] = 1;
			  children_array[2] = 2;
			  break;
		case 8:   n_children = 1;
		          children_array[0] = 3;
			  break;
		case 9:   n_children = 2;
                          children_array[0] = 0;
                          children_array[1] = 3;
			  break;				    
		case 10:  n_children = 2;
                          children_array[0] = 1;
			  children_array[1] = 3;
			  break;
                case 11:  n_children = 3;
		          children_array[0] = 0;
			  children_array[1] = 1;
			  children_array[2] = 3;
			  break;
		case 12:  n_children = 2;
		          children_array[0] = 2;
			  children_array[1] = 3;
			  break;
		case 13:  n_children = 3;
		          children_array[0] = 0;
			  children_array[1] = 2;
			  children_array[2] = 3;
			  break;
		case 14:  n_children = 3;
		          children_array[0] = 1;
			  children_array[1] = 2;
			  children_array[2] = 3;
			  break;
		case 15:  n_children = 4;
		          children_array[0] = 0;
			  children_array[1] = 1;
			  children_array[2] = 2;
			  children_array[3] = 3;
                          break;
	    }
            else 
            if (k_d == 2)
            switch(bv[level].get_2_bits(node)) {
                case 0: n_children = 0;
                        break;
                case 1: n_children = 1;
                        children_array[0] = 0;
                        break;
                case 2: n_children = 1;
                        children_array[0] = 1;
                        break;
                case 3: n_children = 2;
                        children_array[0] = 0;
                        children_array[1] = 1;
                        break;
            }

        }	

        
        uint64_t total_ones_level(uint16_t level)
        {
            return total_ones[level];        
        }        
        

        uint64_t getKD()
        {
            return k_d;        
        }
                
        
        uint16_t getHeight()
        {
            return height;        
        }

                 
 /*       void print(std::ofstream &ost) {         
           uint64_t i, j, aa_r, zz;
                      
           for (i = 0; i < height; i++) {
               if (bv[i].size() > 0)
                   aa_r = bv_rank[i].rank(bv[i].size());
               else
                   aa_r = 0;               
               ost << aa_r << ": ";
            	
               for (zz = j = 0; zz < aa_r; j++) {
                   if ((bv[i])[j]) {
                       ost << '1';
                       zz++;                   
                   }
                   else
                       ost << '0';
               }
               ost << endl;         
            }
        }*/

};
#endif
