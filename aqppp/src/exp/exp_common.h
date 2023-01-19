#pragma once
#include <assert.h>

#include <fstream>
#include <vector>

#include "../aqpp/common_content.h"
#include "../aqpp/sql_interface.h"
#include "../aqpp/tool.h"

namespace expDemo {
void WritePar(FILE* par_file, aqppp::Settings par);

std::pair<double, double> ReadSamples(
    SQLHANDLE& sqlconnectionhandle, aqppp::Settings& o_PAR, int run_num,
    std::vector<std::vector<double>>& o_sample,
    std::vector<std::vector<double>>& o_small_sample);

// pair<double, double> query_real_value4sum(vector<Condition>& cur_q, string
// table_name, MYSQL* conn, const Settings PAR);

std::string FormQueryString(
    const std::vector<aqppp::Condition>& conditions,
    const std::vector<std::string> condition_column_names,
    const std::string aggregation, const std::string aggregate_column_name,
    const std::string db_name, const std::string table_name,
    std::vector<std::unordered_map<int, std::string>>& distinct_itos);
const int QueryRealValue(const std::vector<aqppp::Condition>& conditions,
                         const std::vector<std::string> condition_column_names,
                         const std::string aggregation,
                         const std::string aggregate_column_name,
                         const std::string db_name,
                         const std::string table_name,
                         SQLHANDLE& sql_connection_handle, double& value,
                         bool CLEANCACHE = false);

// double count_db(string table_name, SQLHANDLE &sqlstatementhandle, string
// db_name);

bool IsFileExists(const std::string& filename);

void ReadDirectQueries(std::string direct_file_name, const int DIM,
                       std::vector<std::vector<double>>& o_direct_query_res);

std::vector<double> GetDirectQueries(
    int qid, const std::vector<aqppp::Condition>& cur_query,
    const std::vector<std::vector<double>>& direct_query_res);

static void fprintSample(FILE* out_file,
                         std::vector<std::vector<double>> sample,
                         std::string start_inf, std::string end_inf,
                         int row_limit = INT_MAX);

}  // namespace expDemo