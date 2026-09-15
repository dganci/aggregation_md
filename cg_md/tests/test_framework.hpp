#pragma once


#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace cg::test {

struct Case {
    std::string name;
    std::function<void()> fn;
};

inline std::vector<Case>& registry() {
    static std::vector<Case> cases;
    return cases;
}

struct Registrar {
    Registrar(std::string name, std::function<void()> fn) {
        registry().push_back({std::move(name), std::move(fn)});
    }
};

struct AssertionFailure {
    std::string message;
};

} // namespace cg::test

#define CG_TEST_CONCAT_(a, b) a##b
#define CG_TEST_CONCAT(a, b) CG_TEST_CONCAT_(a, b)

#define CG_TEST(name)                                                                 \
    static void CG_TEST_CONCAT(cg_test_fn_, name)();                                  \
    static const ::cg::test::Registrar CG_TEST_CONCAT(cg_test_registrar_, name)(       \
        #name, CG_TEST_CONCAT(cg_test_fn_, name));                                    \
    static void CG_TEST_CONCAT(cg_test_fn_, name)()

#define CG_CHECK(cond)                                                                \
    do {                                                                              \
        if (!(cond)) {                                                                \
            throw ::cg::test::AssertionFailure{                                       \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) +              \
                ": CG_CHECK failed: " #cond};                                         \
        }                                                                             \
    } while (0)

#define CG_CHECK_EQ(a, b)                                                             \
    do {                                                                              \
        if (!((a) == (b))) {                                                          \
            throw ::cg::test::AssertionFailure{                                       \
                std::string(__FILE__) + ":" + std::to_string(__LINE__) +              \
                ": CG_CHECK_EQ failed: " #a " == " #b};                               \
        }                                                                             \
    } while (0)
