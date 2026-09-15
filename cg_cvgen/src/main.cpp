#include "Config.hpp"
#include "Runner.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    try {
        return cgcv::Runner(cgcv::parse_args(argc, argv)).run();
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        std::cerr << "Use --help for usage.\n";
        return 1;
    }
}
