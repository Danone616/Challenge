#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<math.h>
#include<chrono>
#include<thread>
#include<mutex>

using namespace std;
using namespace chrono;
const double BLOCK_SIZE = 0.051;
const int NUM_THREADS = 16;
struct point {
    double x, y, z;
};
double dist(const point& p1, const point& p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}
double xmi = 1e9, ymi = 1e9, zmi = 1e9, xma = -1e9, yma = -1e9, zma = -1e9, xs, ys, zs;
double a, b, c;
int n, nx, ny, nz;

inline int pos(int x, int y, int z) {
    return x * ny * nz + y * nz + z;
}
vector<point> points;
int res = 0, res1 = 0;
mutex mut;
void func(int istart, int ifinish, vector<vector<int>*>& mat) {
    if (ifinish > nx + 1)ifinish = nx - 1;
    if (istart == 0)istart++;

    for (int ix = istart; ix < ifinish; ix++) {
        for (int jy = 1; jy < ny - 1; jy++) {
            for (int kz = 1; kz < nz - 1; kz++) {

                if (mat[pos(ix, jy, kz)] == NULL)continue;
                int cn = mat[pos(ix, jy, kz)]->size();

                for (int i = 0; i < cn; i++) {
                    for (int j = i + 1; j < cn; j++) {
                        if (dist(points[(*mat[pos(ix, jy, kz)])[i]], points[(*mat[pos(ix, jy, kz)])[j]]) <= 0.05) {
                            lock_guard<mutex> lock(mut);
                            res++;
                        }
                    }
                }
                //k=0
                for (int i = 0; i <= 1; i++) {
                    for (int j = -1; j <= 1; j++) {
                        if (i == 0 && j == -1)j += 2;
                        if (mat[pos(ix + i, jy + j, kz)] == NULL)continue;
                        for (auto& p : *mat[pos(ix, jy, kz)]) {
                            for (auto& q : *mat[pos(ix + i, jy + j, kz)]) {
                                if (dist(points[p], points[q]) <= 0.05) {
                                    lock_guard<mutex> lock(mut);
                                    res1++;
                                }
                            }
                        }
                    }
                }
                //k=1
                for (int i = -1; i <= 1; i++) {
                    for (int j = -1; j <= 1; j++) {
                        if (mat[pos(ix + i, jy + j, kz + 1)] == NULL)continue;
                        for (auto& p : *mat[pos(ix, jy, kz)]) {
                            for (auto& q : *mat[pos(ix + i, jy + j, kz + 1)]) {
                                if (dist(points[p], points[q]) <= 0.05) {
                                    lock_guard<mutex> lock(mut);
                                    res1++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int main() {
    cout << nx << ' ' << ny << ' ' << nz << endl; auto ensq = high_resolution_clock::now();
    ifstream in("In.txt");
    while (in >> a >> b >> c) {
        point p;
        p.x = a;
        p.y = b;
        p.z = c;
        points.push_back(p);
        xmi = min(xmi, a);
        ymi = min(ymi, b);
        zmi = min(zmi, c);
        xma = max(xma, a);
        yma = max(yma, b);
        zma = max(zma, c);
    }
    auto start = high_resolution_clock::now();
    n = points.size();
    cout << n << endl;
    xs = xma - xmi;
    ys = yma - ymi;
    zs = zma - zmi;
    cout << xs << ' ' << ys << ' ' << zs << endl;

    nx = xs / BLOCK_SIZE + 2;
    ny = ys / BLOCK_SIZE + 2;
    nz = zs / BLOCK_SIZE + 2;

    auto read_finish = high_resolution_clock::now();
    vector<vector<int>*>mat(nx * ny * nz, NULL);

    for (int i = 0; i < n; i++) {
        int ix = (points[i].x - xmi) / BLOCK_SIZE + 1;
        int iy = (points[i].y - ymi) / BLOCK_SIZE + 1;
        int iz = (points[i].z - zmi) / BLOCK_SIZE + 1;
        if (mat[pos(ix, iy, iz)] == NULL)mat[pos(ix, iy, iz)] = new vector<int>;
        mat[pos(ix, iy, iz)]->push_back(i);
    }

    auto memory_init = high_resolution_clock::now();

    vector<thread> thread_pool;
    int size_of_block = (nx - 1) / NUM_THREADS + 1;
    for (int i = 0; i < NUM_THREADS; i++) {
        int istart = i * size_of_block;
        int ifinish = (i + 1) * size_of_block;
        thread_pool.push_back(thread(func, istart, ifinish, ref(mat)));
    }
    for (int i = 0; i < thread_pool.size(); i++)thread_pool[i].join();
    auto end = high_resolution_clock::now();
    cout << res << ' ' << res1 << ' ' << res + res1 << endl;
    auto dur_read = read_finish - start;
    auto dur_memory = memory_init - read_finish;
    auto dur_algo = end - memory_init;
    auto dur_total = end - start;

    cout << "read: " << dur_read.count() / 1000000000.0 << ' ' << "sec" << endl;
    cout << "memory init: " << dur_memory.count() / 1000000000.0 << ' ' << "sec" << endl;
    cout << "algo: " << dur_algo.count() / 1000000000.0 << ' ' << "sec" << endl;
    cout << "total: " << dur_total.count() / 1000000000.0 << ' ' << "sec" << endl;

}