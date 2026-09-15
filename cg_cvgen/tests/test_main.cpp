#include "test_framework.hpp"

#include <exception>
#include <iostream>

int main() {
    int failed = 0;
    for (const auto& c : cgcv::test::registry()) {
        try {
            c.fn();
            std::cout << "[PASS] " << c.name << '\n';
        } catch (const cgcv::test::AssertionFailure& e) {
            std::cout << "[FAIL] " << c.name << ": " << e.message << '\n';
            ++failed;
        } catch (const std::exception& e) {
            std::cout << "[FAIL] " << c.name << ": unexpected exception: " << e.what() << '\n';
            ++failed;
        } catch (...) {
            std::cout << "[FAIL] " << c.name << ": unexpected non-standard exception\n";
            ++failed;
        }
    }

    const auto total = cgcv::test::registry().size();
    std::cout << '\n' << (total - static_cast<std::size_t>(failed)) << "/" << total << " tests passed\n";
    return failed == 0 ? 0 : 1;
}
