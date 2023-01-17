// This should contain an entry-point function or class called "RunComparison"
// that is similar to the RunOnce function in comprehensixe_exp, but just
// evaluates the queries that are relevant for my experiments. There should
// probably also be a function that sets the parameters.

#include "comparison.h"

namespace exp_comparison {

const std::vector<std::string> ComparisonExperiment::AGGREGATIONS = {
    "COUNT", "SUM", "AVG"};

const ComparisonExperiment::Parameters ComparisonExperiment::LoadParameters() {
  // TODO load parameter values from text file
  Parameters par = Parameters();
  par.OUTPUT_PATH = "exp_result/comparison/" + par.DATASET_NAME + "/" +
                    "sample_rate_" + std::to_string(par.SAMPLE_RATE);
  par.QUERIES_PATH = par.QUERIES_BASE_PATH + "/" + par.DATASET_NAME + "/" +
                     par.QUERIES_FILENAME;
  return par;
}

// Loads all data from a given database table. Returns the duration.
const double ComparisonExperiment::ReadDB(
    SQLHANDLE& sql_connection_handle, std::vector<std::vector<double>>& o_table,
    const std::string db_name, const std::string table_name) {
  // Setup
  double t_start = clock();
  o_table = std::vector<std::vector<double>>();

  // Get table column names
  std::vector<std::string> column_names =
      aqppp::SqlInterface::GetTableColumns(sql_connection_handle, table_name)
          .all_columns;
  short int n_columns = column_names.size();

  // Define query
  SQLHANDLE sqlstatementhandle = NULL;
  if (SQLAllocHandle(SQL_HANDLE_STMT, sql_connection_handle,
                     &sqlstatementhandle) != SQL_SUCCESS)
    return -1;
  std::string query = "SELECT * FROM " + db_name + "." + table_name + ";";
  aqppp::SqlInterface::SqlQuery(query, sqlstatementhandle);

  // Initialise the output vector by adding a vector for each database column
  // SQLNumResultCols(sqlstatementhandle, &n_columns);
  for (int i = 0; i < n_columns; i++) o_table.push_back(std::vector<double>());

  // Loop over database until all data is added to the return vector
  while (SQLFetch(sqlstatementhandle) == SQL_SUCCESS) {
    for (int i = 0; i < n_columns; i++) {
      double data = aqppp::SqlInterface::Column2Numeric(sqlstatementhandle, i,
                                                        column_names[i]);
      o_table[i].push_back(data);
    }
  }
  SQLFreeHandle(SQL_HANDLE_STMT, sqlstatementhandle);
  double read_time = ((double)clock() - t_start) / CLOCKS_PER_SEC;
  return read_time;
}

// Load queries
const int ComparisonExperiment::LoadQueries(
    const std::string query_filepath, std::vector<aqppp::Query>& o_queries,
    const int n_columns) {
  std::ifstream query_file(query_filepath);
  std::string line;
  o_queries = std::vector<aqppp::Query>();  // init empty
  if (query_file.is_open()) {
    std::getline(query_file, line);  // skip first line (headers)
    while (std::getline(query_file, line)) {
      // For each row, get the condition column id, low bound and upper bound
      std::stringstream ss(line);
      std::string str;
      std::getline(ss, str, ',');
      int condition_column_id = std::stoi(str);
      std::getline(ss, str, ',');
      double bound_low = std::stod(str);
      std::getline(ss, str, ',');
      double bound_high = std::stod(str);

      // Create a query for each combination of aggregation column and
      // aggregation
      for (int agg_col_id = 0; agg_col_id < n_columns; agg_col_id++) {
        for (auto aggregation : AGGREGATIONS) {
          aqppp::Condition condition;
          condition.lb = bound_low;
          condition.ub = bound_high;
          aqppp::Query query = {
              agg_col_id, aggregation, {condition_column_id}, {condition}};
          o_queries.push_back(query);
        }
      }
    }
  } else {
    return -1;
  }
  query_file.close();
  return 0;
}

// Returns the median time to retrieve the sample and small sample database
// tables. The medians are computed over n_runs iterations.
const std::pair<double, double> ComparisonExperiment::ReadSamples(
    SQLHANDLE& sql_connection_handle, const std::string db_name,
    const std::string sample_table_name,
    const std::string small_sample_table_name,
    std::vector<std::vector<double>>& o_sample,
    std::vector<std::vector<double>>& o_small_sample, const int n_runs) {
  std::vector<double> read_sample_times = std::vector<double>();
  std::vector<double> read_small_sample_times = std::vector<double>();
  for (int i = 0; i < n_runs; i++) {
    double time_sample =
        ReadDB(sql_connection_handle, o_sample, db_name, sample_table_name);
    double time_small_sample = ReadDB(sql_connection_handle, o_small_sample,
                                      db_name, small_sample_table_name);
    read_sample_times.push_back(time_sample);
    read_small_sample_times.push_back(time_small_sample);
  }
  double time_read_sample = aqppp::Tool::get_percentile(read_sample_times, 0.5);
  double time_read_small_sample =
      aqppp::Tool::get_percentile(read_small_sample_times, 0.5);
  return {time_read_sample, time_read_small_sample};
};

void ComparisonExperiment::WriteParameters(
    FILE* fp, ComparisonExperiment::Parameters par) {
  fprintf(fp, "---- Databases & files -----\n");
  fprintf(fp, "DATASET_NAME        %s\n", par.DATASET_NAME.c_str());
  fprintf(fp, "DB_NAME             %s\n", par.DB_NAME.c_str());
  fprintf(fp, "TABLE_NAME          %s\n", par.TABLE_NAME.c_str());
  fprintf(fp, "SAMPLE_NAME         %s\n", par.SAMPLE_NAME.c_str());
  fprintf(fp, "SUB_SAMPLE_NAME     %s\n", par.SUB_SAMPLE_NAME.c_str());
  fprintf(fp, "QUERIES_FILENAME    %s\n", par.QUERIES_FILENAME.c_str());
  fprintf(fp, "OUTPUT_PATH         %s\n", par.OUTPUT_PATH.c_str());
  fprintf(fp, "QUERIES_PATH        %s\n", par.QUERIES_PATH.c_str());
  fprintf(fp, "----------------------------\n\n");

  fprintf(fp, "---- Sampling & queries ----\n");
  fprintf(fp, "SAMPLE_RATE         %.3f\n", par.SAMPLE_RATE);
  fprintf(fp, "SUB_SAMPLE_RATE     %.3f\n", par.SUB_SAMPLE_RATE);
  fprintf(fp, "RAND_SEED           %.3f\n", par.RAND_SEED);
  fprintf(fp, "CI_INDEX            %.3f\n", par.CI_INDEX);
  fprintf(fp, "SAMPLE_ROW_NUM      %.3f\n", par.SAMPLE_ROW_NUM);
  fprintf(fp, "NF_MAX_ITER         %.3f\n", par.NF_MAX_ITER);
  fprintf(fp, "INIT_DISTINCT_EVEN  %i\n", par.INIT_DISTINCT_EVEN);
  fprintf(fp, "ALL_MTL_POINTS      %i\n", par.ALL_MTL_POINTS);
  fprintf(fp, "EP_PIECE_NUM        %i\n", par.EP_PIECE_NUM);
  fprintf(fp, "----------------------------\n");
}

// Output files include:
// - parameters (same as loaded from input parameters file)
// - info (timings and high-level statistics)
// - results (query results)
const int ComparisonExperiment::RunExperiment(
    SQLHANDLE& sql_connection_handle) {
  // Setup
  double t_start = clock();
  ComparisonExperiment::Parameters PAR = LoadParameters();
  FILE *parameters_file, *info_file, *results_file, *log_file;
  aqppp::Tool::MkDirRecursively(PAR.OUTPUT_PATH);  // ensure output path exists
  fopen_s(&parameters_file, (PAR.OUTPUT_PATH + "/parameters.txt").data(), "w");
  fopen_s(&info_file, (PAR.OUTPUT_PATH + "/info.txt").data(), "w");
  fopen_s(&results_file, (PAR.OUTPUT_PATH + "/results.txt").data(), "w");
  fopen_s(&log_file, (PAR.OUTPUT_PATH + "/log.txt").data(), "w");

  // Export parameters
  ComparisonExperiment::WriteParameters(parameters_file, PAR);
  fclose(parameters_file);

  // Generate sample database tables and load into memory
  std::pair<double, double> time_samples_creation =
      aqppp::SqlInterface::CreateDbSamples(
          sql_connection_handle, PAR.RAND_SEED, PAR.DB_NAME, PAR.TABLE_NAME,
          {PAR.SAMPLE_RATE, PAR.SUB_SAMPLE_RATE},
          {PAR.SAMPLE_NAME, PAR.SUB_SAMPLE_NAME});
  std::vector<std::vector<double>> sample;
  std::vector<std::vector<double>> small_sample;
  std::pair<double, double> time_read_samples =
      ReadSamples(sql_connection_handle, PAR.DB_NAME, PAR.SAMPLE_NAME,
                  PAR.SUB_SAMPLE_NAME, sample, small_sample);
  PAR.SAMPLE_ROW_NUM = sample[0].size();
  PAR.N_COLUMNS = sample.size();

  // Load queries
  std::vector<aqppp::Query> queries;
  ComparisonExperiment::LoadQueries(PAR.QUERIES_PATH, queries, PAR.N_COLUMNS);

  // Generate prefix cube
  // TODO currently this only works for SUM aggregations. In fact, the "sum"
  // argument in GetPrefixSumCube doesn't actually do anything. Maybe I can
  // extend that function to support other kinds of aggregations. It should
  // be relatively simple if the function is already filtering the data for
  // each cube. I would just need to add more aggregations functions.
  // Possibly store a list of aggregation functions in settings/parameters.
  // Note that more prefix cubes increase the storage requirement.
  //
  // how does it decide how large a sample to make and how big a prefix cube
  // to make? Is there a given space budget at the start? Then it uses a
  // hill-climbing optimisation approach to select the best partition scheme
  // Outcome: prefix cube

  // CAsample is diferent for each aggregation column

  // iterate over all columns (aggregation columns)
  // compute CA sample (temporary, only needed while creating the prefix cube)
  // assign budget using the CAsample and mtl_nums (also only used inside this
  // loop) hill climbing to choose NF_mtl_points (where to set the boundaries in
  // the cube) precompute

  // thus, maybe I don't need d^2 prefix cubes and instead only d
  // this will work if I can just omit conditions on all but one dimension
  // when it comes to queries

  // Compute prefix cubes for AQP. A unique cube is computed for each possible
  // aggregation column. The resulting cubes and resolutions are stored in
  // vectors defined below, which are later accessed for running queries.
  double time_transform_sample = 0;
  double time_prefix_cube_prepare = 0;
  double time_prefix_cube_compute = 0;
  std::vector<std::vector<std::vector<aqppp::CA>>> NF_mtl_points_all;
  std::vector<aqppp::MTL_STRU> NF_mtl_res_all;
  std::vector<std::string> column_names =
      aqppp::SqlInterface::GetTableColumns(sql_connection_handle,
                                           PAR.TABLE_NAME)
          .all_columns;
  for (int agg_col_id = 0; agg_col_id < PAR.N_COLUMNS; agg_col_id++) {
    // Compute CAsample (only needed for computing the prefix cube)
    double t_trans_start = clock();
    std::vector<std::vector<aqppp::CA>> CAsample =
        std::vector<std::vector<aqppp::CA>>();
    aqppp::Tool::TransSample(sample, CAsample, agg_col_id);
    time_transform_sample += (clock() - t_trans_start) / CLOCKS_PER_SEC;

    // Allocate space budget for prefix cube
    double t_prepare_start = clock();
    std::vector<int> mtl_nums;
    aqppp::AssignBudgetForDimensions(PAR.SAMPLE_RATE, PAR.ALL_MTL_POINTS,
                                     PAR.EP_PIECE_NUM, PAR.SAMPLE_ROW_NUM,
                                     PAR.CI_INDEX, PAR.NF_MAX_ITER,
                                     PAR.INIT_DISTINCT_EVEN, false)
        .AssignBudget(CAsample, mtl_nums);

    // Select the critical points for the prefix cube
    std::vector<std::vector<aqppp::CA>> NF_mtl_points =
        std::vector<std::vector<aqppp::CA>>();
    std::vector<double> max_errs;
    std::vector<int> iter_nums;
    aqppp::HillClimbing(PAR.SAMPLE_ROW_NUM, PAR.SAMPLE_RATE, PAR.CI_INDEX,
                        PAR.NF_MAX_ITER, PAR.INIT_DISTINCT_EVEN)
        .ChoosePoints(CAsample, mtl_nums, NF_mtl_points, max_errs, iter_nums);
    time_prefix_cube_prepare += (clock() - t_prepare_start) / CLOCKS_PER_SEC;

    // Compute prefix cube
    std::string aggregate_column = column_names[agg_col_id];
    std::unordered_set<std::string> condition_columns(column_names.begin(),
                                                      column_names.end());
    condition_columns.erase(aggregate_column);
    aqppp::MTL_STRU NF_mtl_res = aqppp::MTL_STRU();
    // This is not used in the "comprehensive experiment".
    if (PAR.isMTL) {
      time_prefix_cube_compute +=
          aqppp::Precompute(PAR.DB_NAME, PAR.TABLE_NAME, aggregate_column,
                            condition_columns)
              .GetPrefixSumCube(NF_mtl_points, sql_connection_handle,
                                NF_mtl_res, "sum");
    }

    // Record prefix cube
    NF_mtl_res_all.push_back(NF_mtl_res);
    NF_mtl_points_all.push_back(NF_mtl_points);
  }

  // Evaluate queries
  // TODO support additional aggregations
  double t_before_queries = clock();
  double total_time_sampling = 0;
  double total_time_aqppp = 0;
  double total_time_exact = 0;
  double n_queries_completed = 0;
  std::vector<double> errors_sampling = std::vector<double>();
  std::vector<double> errors_aqppp = std::vector<double>();
  // TODO add timing to results file
  fprintf(
      results_file,
      "query_id,estimate_sampling,estimate_aqppp,ci_sampling,ci_aqppp,exact_"
      "value,time_sampling,time_aqppp,time_exact,error_sampling,error_aqppp\n");
  std::cout << "Evaluating queries..." << std::endl;
  for (int query_id = 0; query_id < queries.size(); query_id++) {
    // Sampling only
    double t_sampling_start = clock();
    std::pair<double, double> result_sampling = {0, 0};
    // std::pair<double, double> result_sampling =
    //     aqppp::Sampling(PAR.SAMPLE_RATE, PAR.CI_INDEX)
    //         .SamplingForSumQuery(sample, queries[query_id]);
    double duration_sampling = (clock() - t_sampling_start) / CLOCKS_PER_SEC;

    // AQP++
    double t_aqppp_start = clock();
    std::pair<double, double> result_aqppp = {0, 0};
    // std::pair<double, double> result_aqppp =
    //     aqppp::Aqpp(PAR.SAMPLE_ROW_NUM, PAR.SAMPLE_RATE, PAR.CI_INDEX)
    //         .AqppSumQuery(query_id, info_file, sample, small_sample,
    //                       NF_mtl_points, NF_mtl_res, queries[query_id]);
    double duration_aqppp = (clock() - t_aqppp_start) / CLOCKS_PER_SEC;

    // Exact value... not necessary since I compute this elsewhere already?
    double t_exact_start = clock();
    double exact_value = 0;
    // double exact_value = expDemo::QueryRealValue(
    //     queries[query_id], PAR.TABLE_NAME, sql_connection_handle, PAR,
    //     "sum");
    double duration_exact = (clock() - t_exact_start) / CLOCKS_PER_SEC;

    // Compute errors and store results
    // TODO check if these are actually errors... I think it is just the ratio
    // between the confidence interval width and the exact value
    double error_sampling =
        exact_value > 0 ? result_sampling.second / exact_value : 0;
    double error_aqppp =
        exact_value > 0 ? result_aqppp.second / exact_value : 0;
    errors_sampling.push_back(error_sampling);
    errors_aqppp.push_back(error_aqppp);
    total_time_sampling += duration_sampling;
    total_time_aqppp += duration_aqppp;
    total_time_exact += duration_exact;
    fprintf(results_file, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", query_id,
            result_sampling.first, result_aqppp.first, result_sampling.second,
            result_aqppp.second, exact_value, duration_sampling, duration_aqppp,
            duration_exact, error_sampling, error_aqppp);

    n_queries_completed++;
  }
  fclose(results_file);
  fclose(log_file);

  // Export everything
  double error_sampling_mean = aqppp::Tool::get_avg(errors_sampling);
  double error_aqppp_mean = aqppp::Tool::get_avg(errors_aqppp);
  double error_sampling_median =
      aqppp::Tool::get_percentile(errors_sampling, 0.5);
  double error_aqppp_median = aqppp::Tool::get_percentile(errors_aqppp, 0.5);
  double error_sampling_99p =
      aqppp::Tool::get_percentile(errors_sampling, 0.99);
  double error_aqppp_99p = aqppp::Tool::get_percentile(errors_aqppp, 0.99);
  double time_preprocessing = time_samples_creation.first +
                              time_samples_creation.second +
                              time_transform_sample + time_prefix_cube_prepare +
                              time_prefix_cube_compute;
  double time_before_queries = (t_before_queries - t_start) / CLOCKS_PER_SEC;
  double time_read_sample = time_read_samples.first;
  double time_read_small_sample = time_read_samples.second;
  double time_exp = (clock() - t_start) / CLOCKS_PER_SEC;

  fprintf(info_file, "Total experiment duration    %.3f s\n", time_exp);
  fprintf(info_file, "Sample generation time       %.3f s\n",
          time_samples_creation.first);
  fprintf(info_file, "Subsample generation time    %.3f s\n",
          time_samples_creation.second);
  fprintf(info_file, "\n");

  fprintf(info_file, "Total sample rows            %i\n", PAR.SAMPLE_ROW_NUM);
  fprintf(info_file, "\n");

  fprintf(info_file, "Total pre-processing time    %.3f s\n",
          time_preprocessing);
  fprintf(info_file, "Total pre-query time         %.3f s\n",
          time_before_queries);
  fprintf(info_file, "Time to read sample          %.3f s\n", time_read_sample);
  fprintf(info_file, "Time to read sub-sample      %.3f s\n",
          time_read_small_sample);
  fprintf(info_file, "Average query time sampling  %.3f s\n",
          total_time_sampling / n_queries_completed + time_read_sample);
  fprintf(info_file, "Average query time AQP++     %.3f s\n",
          total_time_aqppp / n_queries_completed + time_read_sample +
              time_read_small_sample);
  fprintf(info_file, "Average query time exact     %.3f s\n",
          total_time_exact / n_queries_completed);
  fprintf(info_file, "\n");

  fprintf(info_file, "Average error sampling 	     %.3f %%\n",
          100 * error_sampling_mean);
  fprintf(info_file, "Average error AQP++          %.3f %%\n",
          100 * error_aqppp_mean);
  fprintf(info_file, "Median error sampling 	     %.3f %%\n",
          100 * error_sampling_median);
  fprintf(info_file, "Median error AQP++ 	         %.3f %%\n",
          100 * error_aqppp_median);
  fprintf(info_file, "99th percentile sampling     %.3f %%\n",
          100 * error_sampling_99p);
  fprintf(info_file, "99th percentile AQP++ 	     %.3f %%\n",
          100 * error_aqppp_99p);
  fprintf(info_file, "\n");

  // NOTE: according to the paper, they present only SUM queries.
  // They claim that it can be adapted for COUNT and AVERAGE queries easily.
  // (Possibly take the same approach as DeepDB by computing the sum and count
  // queries and using these to compute the average. I think DeepDB did count
  // and average and then derived sum.) Actually, they do have some
  // suggestions in the paper, including creating a summy column for count
  // queries... maybe this can be done better... They note that creating a
  // prefix cube for these is a bit more complex and give an explanation for
  // each... which I probably need to implement.

  fclose(info_file);
  std::cout << "Comparison experiment completed." << std::endl;
  return 0;
}
}  // namespace exp_comparison