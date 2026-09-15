#include "test_framework.hpp"
#include "Config.hpp"
#include "FileUtils.hpp"

#include <string>

using namespace cg;

namespace {

Config base() {
    Config c;
    c.project_dir = "/tmp/pj";
    c.n_prot = 5;
    c.protomer_name = "desmin_head";
    c.seq_length = 108;
    c.atoms_per_prot = 237;
    c.use_elastic = false;
    c.seed = 20260729;
    c.md_dt_ps = 0.01;
    c.adaptive = true;
    c.adaptive_chunk_us = 0.0005;
    c.adaptive_max_total_us = 0.001;
    c.adaptive_min_total_us = 0.0005;
    return c;
}

} // namespace

CG_TEST(identical_parameters_give_the_same_run_id) {
    CG_CHECK_EQ(base().runId(), base().runId());
    CG_CHECK_EQ(base().runName(), std::string("5xdesmin_head_") + base().runId());
}

CG_TEST(physics_changes_fork_a_new_run) {
    const auto ref = base().runId();

    auto hotter = base();          hotter.temperature_K = 320.0;
    auto saltier = base();         saltier.salt_M = 0.30;
    auto bigger = base();          bigger.n_prot = 10;
    auto other_pdb = base();       other_pdb.protomer_name = "desmin_1-263";
    auto phos = base();            phos.is_phospho = true;
    auto reseeded = base();        reseeded.seed = 1;
    auto elastic = base();         elastic.use_elastic = true;
    auto ss = base();              ss.ss_string = std::string(108, 'H');
    auto dssp = base();            dssp.use_dssp = true;
    auto relaxed = base();         relaxed.relax = true;
    auto eps = base();             eps.epsilon_r = 12.0;
    auto thermo = base();          thermo.thermostat_mode = "legacy";
    auto stride = base();          stride.plumed_stride = 500;

    for (const auto& variant : {hotter, saltier, bigger, other_pdb, phos, reseeded,
                                elastic, ss, dssp, relaxed, eps, thermo, stride}) {
        CG_CHECK(variant.runId() != ref);
    }
}

CG_TEST(segment_length_forks_a_new_run) {
    auto shorter = base();
    shorter.adaptive_chunk_us = 0.0001;
    CG_CHECK(shorter.runId() != base().runId());
}

CG_TEST(execution_knobs_do_not_fork_a_run) {
    const auto ref = base().runId();

    auto threads = base();   threads.ntomp = 4;
    auto stage = base();     stage.stage = "adaptive";
    auto binary = base();    binary.gmx = "gmx";
    auto pinned = base();    pinned.pin = "on";
    auto quiet = base();     quiet.report = false;
    auto dry = base();       dry.dry_run = true;
    auto elsewhere = base(); elsewhere.project_dir = "/other";

    for (const auto& variant : {threads, stage, binary, pinned, quiet, dry, elsewhere}) {
        CG_CHECK_EQ(variant.runId(), ref);
    }
}

CG_TEST(length_and_stop_rules_do_not_fork_a_run) {
    const auto ref = base().runId();

    auto longer = base();       longer.adaptive_max_total_us = 0.005;
    auto sooner = base();       sooner.adaptive_min_total_us = 0.0005;
    auto looser = base();       looser.adaptive_min_cn_range = 99.0;
    auto patterns = base();     patterns.adaptive_min_unique_contact_patterns = 99;
    auto md_total = base();     md_total.md_total_us = 42.0;

    for (const auto& variant : {longer, sooner, looser, patterns, md_total}) {
        CG_CHECK_EQ(variant.runId(), ref);
    }
}

CG_TEST(run_tag_overrides_the_derived_id) {
    auto tagged = base();
    tagged.run_tag = "pilot_A";
    CG_CHECK_EQ(tagged.runId(), std::string("pilot_A"));
    CG_CHECK_EQ(tagged.runName(), std::string("5xdesmin_head_pilot_A"));
    CG_CHECK_EQ(tagged.runDir(), tagged.project_dir / "runs" / "5xdesmin_head_pilot_A");
}

CG_TEST(all_prep_artifacts_live_under_the_run_directory) {
    const auto c = base();
    const auto run = c.runDir().string();
    for (const auto& p : {c.cgPath(), c.cgRelaxedPath(), c.topologyPath(), c.packRawPath(),
                          c.packCleanPath(), c.packInputPath(), c.boxedPath(),
                          c.solvatedPath(), c.beadsPath(), c.systemDir(), c.resultDir()}) {
        CG_CHECK(p.string().rfind(run, 0) == 0);
    }
    CG_CHECK(c.pdbPath().string().rfind(run, 0) != 0);
}

CG_TEST(identity_string_is_readable_and_covers_the_key_fields) {
    const auto text = base().identityString();
    for (const auto* key : {"n_prot = 5", "protomer = desmin_head", "temperature_K",
                            "adaptive_chunk_us", "seed = 20260729"}) {
        CG_CHECK(text.find(key) != std::string::npos);
    }
    for (const auto* key : {"ntomp", "stage", "project_dir", "report_script"}) {
        CG_CHECK(text.find(key) == std::string::npos);
    }
}


CG_TEST(lambda_pw_forks_the_run) {
    auto a = base();
    auto b = base();
    b.martini_lambda_pw = 1.10;
    CG_CHECK(a.runId() != b.runId());
}

CG_TEST(lambda_pw_selects_the_rescaled_itp) {
    auto c = base();
    CG_CHECK_EQ(c.martiniMainItp(), std::string("martini_v3.0.0.itp"));
    c.martini_lambda_pw = 1.10;
    CG_CHECK_EQ(c.martiniMainItp(), std::string("martini_v3.0.0_lpw110.itp"));
    c.martini_lambda_pw = 1.06;
    CG_CHECK_EQ(c.martiniMainItp(), std::string("martini_v3.0.0_lpw106.itp"));
}

CG_TEST(lambda_pw_out_of_range_is_refused) {
    auto c = base();
    c.martini_lambda_pw = 3.0;
    bool threw = false;
    try { c.validate(); } catch (const std::exception&) { threw = true; }
    CG_CHECK(threw);
}

CG_TEST(lambda_pw_appears_in_the_identity_string) {
    auto c = base();
    c.martini_lambda_pw = 1.10;
    CG_CHECK(c.identityString().find("martini_lambda_pw = 1.1000") != std::string::npos);
}

CG_TEST(distinct_copies_changes_the_run_identity) {
    Config a = base();
    Config b = base();
    b.distinct_copies = true;
    CG_CHECK(a.runId() != b.runId());
}

CG_TEST(packmol_sources_fall_back_unless_every_copy_exists) {
    const auto root = std::filesystem::temp_directory_path() / "cg_md_test_dcopies";
    std::filesystem::remove_all(root);

    Config c = base();
    c.project_dir = root;
    c.distinct_copies = true;
    std::filesystem::create_directories(c.prepDir());
    write_text(c.cgPath(), "ATOM\n");

    CG_CHECK_EQ(c.packmolSourcePaths().size(), std::size_t{1});
    CG_CHECK_EQ(c.packmolSourcePaths().front(), c.cgPath());

    for (int i = 0; i < c.n_prot - 1; ++i) write_text(c.cgRelaxedCopyPath(i), "ATOM\n");
    CG_CHECK_EQ(c.packmolSourcePaths().size(), std::size_t{1});

    write_text(c.cgRelaxedCopyPath(c.n_prot - 1), "ATOM\n");
    const auto full = c.packmolSourcePaths();
    CG_CHECK_EQ(full.size(), static_cast<std::size_t>(c.n_prot));
    for (int i = 0; i < c.n_prot; ++i) CG_CHECK_EQ(full[i], c.cgRelaxedCopyPath(i));

    c.distinct_copies = false;
    CG_CHECK_EQ(c.packmolSourcePaths().size(), std::size_t{1});

    std::filesystem::remove_all(root);
}

CG_TEST(distinct_copies_without_relax_is_rejected) {
    Config c = base();
    c.distinct_copies = true;
    c.relax = false;
    bool threw = false;
    try { c.validate(); } catch (const std::runtime_error&) { threw = true; }
    CG_CHECK(threw);

    c.relax = true;
    c.relax_us = 0.1;
    bool threw_with_relax = false;
    try { c.validate(); } catch (const std::runtime_error&) { threw_with_relax = true; }
    CG_CHECK(!threw_with_relax);
}
