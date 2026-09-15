#include "test_framework.hpp"
#include "TopologyEditor.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"

#include <exception>
#include <filesystem>

using namespace cg;

namespace {
std::filesystem::path temp_dir(const std::string& name) {
    auto dir = std::filesystem::temp_directory_path() / ("cg_md_test_" + name);
    std::filesystem::create_directories(dir);
    return dir;
}
} // namespace

CG_TEST(find_first_itp_matches_prefix_and_extension) {
    const auto dir = temp_dir("itp_dir");
    write_text(dir / "desmin_head_0.itp", "; stub\n");
    write_text(dir / "unrelated.txt", "not an itp\n");

    CG_CHECK_EQ(find_first_itp(dir, "desmin_head"), "desmin_head_0.itp");

    std::filesystem::remove_all(dir);
}

CG_TEST(find_first_itp_throws_when_nothing_matches) {
    const auto dir = temp_dir("itp_dir_empty");
    bool threw = false;
    try {
        find_first_itp(dir, "nonexistent");
    } catch (const std::exception&) {
        threw = true;
    }
    CG_CHECK(threw);

    std::filesystem::remove_all(dir);
}

CG_TEST(patch_martini_topology_rewrites_includes_and_molecule_count) {
    const auto dir = temp_dir("topo_dir");
    const auto martini_dir = temp_dir("martini_dir");
    const auto topo_path = dir / "protomer_topo.top";

    write_lines(topo_path, {
        "#include \"martini.itp\"",
        "",
        "[ moleculetype ]",
        "; keep this section untouched",
        "",
        "[ molecules ]",
        "; name        count",
        "protomer 1"
    });

    patch_martini_topology(topo_path, martini_dir, "protomer", 10, "protomer_0.itp");

    const auto lines = read_lines(topo_path);

    int include_count = 0;
    bool found_molecule_line = false;
    bool kept_moleculetype_section = false;
    for (const auto& line : lines) {
        if (starts_with(line, "#include")) ++include_count;
        if (line.find("protomer_0") != std::string::npos && line.find("10") != std::string::npos) found_molecule_line = true;
        if (line == "[ moleculetype ]") kept_moleculetype_section = true;
    }

    CG_CHECK_EQ(include_count, 4);
    CG_CHECK(found_molecule_line);
    CG_CHECK(kept_moleculetype_section);

    std::filesystem::remove_all(dir);
    std::filesystem::remove_all(martini_dir);
}

CG_TEST(a_single_molecule_topology_keeps_one_entry_and_drops_the_solvent) {
    const auto dir = std::filesystem::temp_directory_path() / "cg_md_single_topo";
    std::filesystem::create_directories(dir);
    const auto src = dir / "system.top";
    const auto dst = dir / "relax.top";

    write_lines(src, {
        "#include \"martini_v3.0.0.itp\"",
        "",
        "[ system ]",
        "Five protomers in water",
        "",
        "[ molecules ]",
        "; Compound        #mols",
        "protomer_0           5",
        "W                41261",
        "NA+                 89",
        "CL-                 74"
    });

    write_single_molecule_topology(src, dst);
    const auto lines = read_lines(dst);

    int entries = 0;
    bool kept_include = false, kept_comment = false, one_protomer = false;
    bool in_molecules = false;
    for (const auto& line : lines) {
        const auto t = trim(line);
        if (starts_with(t, "#include")) kept_include = true;
        if (t == "; Compound        #mols") kept_comment = true;
        if (starts_with(t, "[")) { in_molecules = starts_with(t, "[ molecules ]"); continue; }
        if (!in_molecules || t.empty() || starts_with(t, ";")) continue;
        ++entries;
        const auto fields = split_ws(t);
        if (fields.size() == 2 && fields[0] == "protomer_0" && fields[1] == "1") one_protomer = true;
    }

    CG_CHECK_EQ(entries, 1);
    CG_CHECK(one_protomer);
    CG_CHECK(kept_include);
    CG_CHECK(kept_comment);
    CG_CHECK_EQ(static_cast<int>(read_lines(src).size()), 11);

    std::filesystem::remove_all(dir);
}
