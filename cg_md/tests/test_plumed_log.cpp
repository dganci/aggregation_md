#include "test_framework.hpp"
#include "PlumedLog.hpp"

#include <string>

using namespace cg;

CG_TEST(walker_count_reads_the_colon_spelling) {
    const std::string log =
        "PLUMED: Action METAD\n"
        "PLUMED:   with label metad\n"
        "PLUMED:   Multiple walkers active using MPI communnication\n"
        "PLUMED:   number of walkers: 4\n"
        "PLUMED:   walker id: 0\n";
    CG_CHECK_EQ(plumed_walker_count(log), 4);
}

CG_TEST(walker_count_reads_the_leading_number_spelling) {
    const std::string log =
        "PLUMED:   Gaussian width  0.100000\n"
        "PLUMED:   4 multiple walkers active\n"
        "PLUMED:   walker id 2\n";
    CG_CHECK_EQ(plumed_walker_count(log), 4);
}

CG_TEST(walker_count_reports_one_when_sharing_is_not_active) {
    CG_CHECK_EQ(plumed_walker_count("PLUMED:   number of walkers: 1\n"), 1);
}

CG_TEST(walker_count_is_unknown_when_the_log_is_silent) {
    CG_CHECK_EQ(plumed_walker_count("PLUMED: Action DISTANCE\nPLUMED:   with label d_1_2\n"), -1);
    CG_CHECK_EQ(plumed_walker_count(""), -1);
}

CG_TEST(walker_count_ignores_prose_that_merely_mentions_walkers) {
    CG_CHECK_EQ(plumed_walker_count("PLUMED:   Multiple walkers active using MPI communnication\n"), -1);
}

CG_TEST(walker_count_takes_the_first_answer_it_finds) {
    const std::string log =
        "PLUMED:   number of walkers: 4\n"
        "PLUMED: ...\n"
        "PLUMED:   number of walkers: 4\n";
    CG_CHECK_EQ(plumed_walker_count(log), 4);
}


CG_TEST(a_numbered_backup_of_the_target_is_recognised) {
    CG_CHECK(is_plumed_backup("bck.0.HILLS", "HILLS"));
    CG_CHECK(is_plumed_backup("bck.12.HILLS", "HILLS"));
    CG_CHECK(is_plumed_backup("bck.last.HILLS", "HILLS"));
    CG_CHECK(is_plumed_backup("bck.3.COLVAR_monitor", "COLVAR_monitor"));
}

CG_TEST(the_live_file_is_not_its_own_backup) {
    CG_CHECK(!is_plumed_backup("HILLS", "HILLS"));
    CG_CHECK(!is_plumed_backup("COLVAR", "COLVAR"));
}

CG_TEST(a_backup_of_a_different_file_does_not_count) {
    CG_CHECK(!is_plumed_backup("bck.0.COLVAR", "HILLS"));
    CG_CHECK(!is_plumed_backup("bck.0.COLVAR_monitor", "COLVAR"));
}

CG_TEST(files_that_merely_look_like_backups_are_ignored) {
    CG_CHECK(!is_plumed_backup("bck.HILLS", "HILLS"));
    CG_CHECK(!is_plumed_backup("bck.zero.HILLS", "HILLS"));
    CG_CHECK(!is_plumed_backup("my_bck.0.HILLS", "HILLS"));
    CG_CHECK(!is_plumed_backup("bck.0.HILLS.old", "HILLS"));
    CG_CHECK(!is_plumed_backup("", "HILLS"));
}
