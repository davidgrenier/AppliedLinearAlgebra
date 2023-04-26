#include <chrono>
#include <ios>

using std::cout;
using std::string;
using namespace std::chrono;

template<typename T>
T bench(string name, const std::function<T()>& f) {
    auto start = high_resolution_clock::now();
    auto t = f();
    auto ellapsed = duration_cast<microseconds>(high_resolution_clock::now()-start).count();

    cout.precision(3);
    cout << name << " took " << std::fixed << ellapsed/1.e3 << std::defaultfloat << "ms\n";

    return t;
}

void bench(string name, const std::function<void()>& f) {
    bench<int>(name, [&] { f(); return 0; });
}
