#include "Config.hpp"
#include "Workflow.hpp"

#include <exception>
#include <iostream>

int main(int argc, char** argv) {
    try {
        cg::Workflow{cg::Config::fromArgs(argc, argv)}.run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << '\n';
        return 1;
    }
}
