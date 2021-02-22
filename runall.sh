cd runqueries/
GR='\033[0;32m'
NC='\033[0m' # No Color
echo -e "Running ${GR}tr1...${NC}"
bash runqueries-tr1-bfs-sorted.sh > outputs/tr1.txt
echo -e "Running ${GR}tr2...${NC}"
bash runqueries-tr2-bfs-sorted.sh > outputs/tr2.txt
echo -e "Running ${GR}j3...${NC}"
bash runqueries-j3-bfs-sorted.sh > outputs/j3.txt
echo -e "Running ${GR}j4...${NC}"
bash runqueries-j4-bfs-sorted.sh > outputs/j4.txt
echo -e "Running ${GR}p2...${NC}"
bash runqueries-p2-bfs-sorted.sh > outputs/p2.txt
echo -e "Running ${GR}p3...${NC}"
bash runqueries-p3-bfs-sorted.sh > outputs/p3.txt
echo -e "Running ${GR}p4...${NC}"
bash runqueries-p4-bfs-sorted.sh > outputs/p4.txt
echo -e "Running ${GR}s1...${NC}"
bash runqueries-s1-bfs-sorted.sh > outputs/s1.txt
echo -e "Running ${GR}s2...${NC}"
bash runqueries-s2-bfs-sorted.sh > outputs/s2.txt
echo -e "Running ${GR}s3...${NC}"
bash runqueries-s3-bfs-sorted.sh > outputs/s3.txt
echo -e "Running ${GR}s4...${NC}"
bash runqueries-s4-bfs-sorted.sh > outputs/s4.txt
echo -e "Running ${GR}t2...${NC}"
bash runqueries-t2-bfs-sorted.sh > outputs/t2.txt
echo -e "Running ${GR}t3...${NC}"
bash runqueries-t3-bfs-sorted.sh > outputs/t3.txt
echo -e "Running ${GR}t4...${NC}"
bash runqueries-t4-bfs-sorted.sh > outputs/t4.txt
echo -e "Running ${GR}ti2...${NC}"
bash runqueries-ti2-bfs-sorted.sh > outputs/ti2.txt
echo -e "Running ${GR}ti3...${NC}"
bash runqueries-ti3-bfs-sorted.sh > outputs/ti3.txt
echo -e "Running ${GR}ti4...${NC}"
bash runqueries-ti4-bfs-sorted.sh > outputs/ti4.txt


