cd queries/
GR='\033[0;32m'
NC='\033[0m' # No Color
echo -e "Building ${GR}tr1...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib tr1.cpp -o ../runqueries/tr1 -lsdsl -pthread
echo -e "Building ${GR}tr2...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib tr2.cpp -o ../runqueries/tr2 -lsdsl -pthread
echo -e "Building ${GR}j3...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib j3.cpp -o ../runqueries/j3 -lsdsl -pthread
echo -e "Building ${GR}j4...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib j4.cpp -o ../runqueries/j4 -lsdsl -pthread
echo -e "Building ${GR}p2...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib p2.cpp -o ../runqueries/p2 -lsdsl -pthread
echo -e "Building ${GR}p3...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib p3.cpp -o ../runqueries/p3 -lsdsl -pthread
echo -e "Building ${GR}p4...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib p4.cpp -o ../runqueries/p4 -lsdsl -pthread
echo -e "Building ${GR}s1...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib s1.cpp -o ../runqueries/s1 -lsdsl -pthread
echo -e "Building ${GR}s2...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib s2.cpp -o ../runqueries/s2 -lsdsl -pthread
echo -e "Building ${GR}s3...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib s3.cpp -o ../runqueries/s3 -lsdsl -pthread
echo -e "Building ${GR}s4...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib s4.cpp -o ../runqueries/s4 -lsdsl -pthread
echo -e "Building ${GR}t2...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib t2.cpp -o ../runqueries/t2 -lsdsl -pthread
echo -e "Building ${GR}t3...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib t3.cpp -o ../runqueries/t3 -lsdsl -pthread
echo -e "Building ${GR}t4...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib t4.cpp -o ../runqueries/t4 -lsdsl -pthread
echo -e "Building ${GR}ti2...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib ti2.cpp -o ../runqueries/ti2 -lsdsl -pthread
echo -e "Building ${GR}ti3...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib ti3.cpp -o ../runqueries/ti3 -lsdsl -pthread
echo -e "Building ${GR}ti4...${NC}"
g++ -std=c++11 -mpopcnt -m64 -frename-registers -O9 -msse4.2 -DNDEBUG -I ~/include -L ~/lib ti4.cpp -o ../runqueries/ti4 -lsdsl -pthread
echo -e "${GR}done building!${NC}"
