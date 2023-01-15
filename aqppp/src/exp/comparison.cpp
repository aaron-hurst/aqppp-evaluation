// This should contain an entry-point function or class called "RunComparison"
// that is similar to the RunOnce function in comprehensixe_exp, but just
// evaluates the queries that are relevant for my experiments. There should
// probably also be a function that sets the parameters.

#include "comparison.h"

namespace exp_comparison {

const std::vector<std::string> AGGREGATIONS = {"COUNT", "SUM", "AVG"};

const ComparisonExperiment::Parameters ComparisonExperiment::LoadParameters() {
  // TODO load parameter values from text file
  Parameters par = Parameters();
  par.OUTPUT_PATH = "exp_result/comparison/" + par.DATASET_NAME + "/" +
                    "sample_rate_" + std::to_string(par.SAMPLE_RATE);
  par.QUERIES_PATH =
      "../../queries/" + par.DATASET_NAME + "/" + par.QUERIES_FILENAME;
  return par;
}

// Loads all data from a given database table. Returns the duration.
const double ComparisonExperiment::ReadDB(
    SQLHANDLE& sql_connection_handle, std::vector<std::vector<double>>& o_table,
    const std::string db_name, const std::string table_name) {
  double t_start = clock();
  o_table = std::vector<std::vector<double>>();

  // Define query
  SQLHANDLE sqlstatementhandle = NULL;
  if (SQLAllocHandle(SQL_HANDLE_STMT, sql_connection_handle,
                     &sqlstatementhandle) != SQL_SUCCESS)
    return -1;
  std::string query = "SELECT * FROM " + db_name + "." + table_name + ";";
  aqppp::SqlInterface::SqlQuery(query, sqlstatementhandle);

  // Initialise the output vector by adding a vector for each database column
  short int n_columns = 0;
  SQLNumResultCols(sqlstatementhandle, &n_columns);
  for (int i = 0; i < n_columns; i++) o_table.push_back(std::vector<double>());

  // Loop over database until all data is added to the return vector
  aqppp::SqlInterface::TableColumns table_columns =
      aqppp::SqlInterface::GetTableColumns(sql_connection_handle, table_name);
  std::vector<std::string> column_names = table_columns.numeric_columns;
  column_names.insert(column_names.end(),
                      table_columns.categorical_columns.begin(),
                      table_columns.categorical_columns.end());
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
    const std::string query_filepath,
    std::vector<ComparisonExperiment::Query> o_queries, const int n_columns) {
  std::ifstream query_file(query_filepath);
  std::string line;
  o_queries = std::vector<ComparisonExperiment::Query>();  // init empty
  while (std::getline(query_file, line)) {
    // For each row, get the condition column id, low bound and upper bound
    std::stringstream ss(line);
    std::string str;
    std::getline(ss, str, ',');
    int condition_column_id = std::stoi(str);
    std::getline(ss, str, ',');
    int bound_low = std::stod(str);
    std::getline(ss, str, ',');
    int bound_high = std::stod(str);

    // Create a query for each combination of aggregation column and aggregation
    for (int agg_col_id = 0; agg_col_id < n_columns; agg_col_id++) {
      for (auto aggregation : AGGREGATIONS) {
        aqppp::Condition condition;
        condition.lb = bound_low;
        condition.ub = bound_high;
        ComparisonExperiment::Query query = {
            agg_col_id, aggregation, {condition_column_id}, {condition}};
        o_queries.push_back(query);
      }
    }
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
  fprintf(fp, "DATASET_NAME       %s\n", par.DATASET_NAME.c_str());
  fprintf(fp, "DB_NAME            %s\n", par.DB_NAME.c_str());
  fprintf(fp, "TABLE_NAME         %s\n", par.TABLE_NAME.c_str());
  fprintf(fp, "SAMPLE_NAME        %s\n", par.SAMPLE_NAME.c_str());
  fprintf(fp, "SUB_SAMPLE_NAME    %s\n", par.SUB_SAMPLE_NAME.c_str());
  fprintf(fp, "QUERIES_FILENAME   %s\n", par.QUERIES_FILENAME.c_str());
  fprintf(fp, "OUTPUT_PATH        %s\n", par.OUTPUT_PATH.c_str());
  fprintf(fp, "QUERIES_PATH       %s\n", par.QUERIES_PATH.c_str());

  fprintf(fp, "---- Sampling & queries ----\n");
  fprintf(fp, "SAMPLE_RATE        %f\n", par.SAMPLE_RATE);
  fprintf(fp, "SUB_SAMPLE_RATE    %f\n", par.SUB_SAMPLE_RATE);
  fprintf(fp, "RAND_SEED          %f\n", par.RAND_SEED);
  fprintf(fp, "CI_INDEX           %f\n", par.CI_INDEX);
  fprintf(fp, "SAMPLE_ROW_NUM     %f\n", par.SAMPLE_ROW_NUM);
  fprintf(fp, "NF_MAX_ITER        %f\n", par.NF_MAX_ITER);
  fprintf(fp, "INIT_DISTINCT_EVEN %f\n", par.INIT_DISTINCT_EVEN);
  fprintf(fp, "ALL_MTL_POINTS     %i\n", par.ALL_MTL_POINTS);
  fprintf(fp, "EP_PIECE_NUM       %i\n", par.EP_PIECE_NUM);
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
  FILE *parameters_file, *info_file, *results_file;
  aqppp::Tool::MkDirRecursively(PAR.OUTPUT_PATH);  // ensure output path exists
  fopen_s(&parameters_file, (PAR.OUTPUT_PATH + "/parameters.txt").data(), "w");
  fopen_s(&info_file, (PAR.OUTPUT_PATH + "/info.txt").data(), "w");
  fopen_s(&results_file, (PAR.OUTPUT_PATH + "/results.txt").data(), "w");

  // Export parameters
  ComparisonExperiment::WriteParameters(parameters_file, PAR);
  fclose(parameters_file);

  // Generate sample database tables, load into memory and generate a CA
  // sample, which is used to generate the prefix cube.
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
  double time_transform_sample = clock();
  std::vector<std::vector<aqppp::CA>> CAsample =
      std::vector<std::vector<aqppp::CA>>();
  aqppp::Tool::TransSample(sample, CAsample);
  time_transform_sample = (clock() - time_transform_sample) / CLOCKS_PER_SEC;

  // Load queries
  std::vector<Query> queries;
  int n_columns = sample[0].size();
  ComparisonExperiment::LoadQueries(PAR.QUERIES_PATH, queries, n_columns);

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
  double t_prefix_cube_start = clock();
  std::vector<std::vector<aqppp::CA>> NF_mtl_points;
  aqppp::MTL_STRU NF_mtl_res = aqppp::MTL_STRU();
  std::vector<int> mtl_nums;
  std::cout << "Start mtl..." << std::endl;
  aqppp::AssignBudgetForDimensions(
      PAR.SAMPLE_RATE, PAR.ALL_MTL_POINTS, PAR.EP_PIECE_NUM, PAR.SAMPLE_ROW_NUM,
      PAR.CI_INDEX, PAR.NF_MAX_ITER, PAR.INIT_DISTINCT_EVEN, false)
      .AssignBudget(CAsample, mtl_nums);
  NF_mtl_points = std::vector<std::vector<aqppp::CA>>();
  std::vector<double> max_errs;
  std::vector<int> iter_nums;
  aqppp::HillClimbing(PAR.SAMPLE_ROW_NUM, PAR.SAMPLE_RATE, PAR.CI_INDEX,
                      PAR.NF_MAX_ITER, PAR.INIT_DISTINCT_EVEN)
      .ChoosePoints(CAsample, mtl_nums, NF_mtl_points, max_errs, iter_nums);
  double time_prefix_cube_prepare =
      (clock() - t_prefix_cube_start) / CLOCKS_PER_SEC;
  double time_prefix_cube_compute =
      aqppp::Precompute(PAR.DB_NAME, PAR.TABLE_NAME, PAR.AGGREGATE_NAME,
                        PAR.CONDITION_NAMES)
          .GetPrefixSumCube(NF_mtl_points, sql_connection_handle, NF_mtl_res,
                            "sum");

  // Evaluate queries
  // TODO support additional aggregations
  double t_before_queries = clock();
  double time_sampling = 0;
  double time_aqppp = 0;
  double time_exact = 0;
  double n_queries_completed = 0;
  std::vector<double> errors_sampling = std::vector<double>();
  std::vector<double> errors_aqppp = std::vector<double>();
  // TODO change format (especially seperator) for results file
  fprintf(results_file,
          "query_id\t spl_est\t spl_ci\t aqpp_est\t aqpp_ci\t real_value\t "
          "time_spl\t time_aqpp\t time_direct\t error_spl\t error_aqpp\n");
  std::cout << "Evaluating queries..." << std::endl;
  for (int query_id = 0; query_id < queries.size(); query_id++) {
    // Sampling only
    double t_sampling_start = clock();
    std::pair<double, double> result_sampling =
        aqppp::Sampling(PAR.SAMPLE_RATE, PAR.CI_INDEX)
            .SamplingForSumQuery(sample, queries[query_id]);
    double duration_sampling = (clock() - t_sampling_start) / CLOCKS_PER_SEC;

    // AQP++
    double t_aqppp_start = clock();
    std::pair<double, double> result_aqppp =
        aqppp::Aqpp(PAR.SAMPLE_ROW_NUM, PAR.SAMPLE_RATE, PAR.CI_INDEX)
            .AqppSumQuery(query_id, info_file, sample, small_sample,
                          NF_mtl_points, NF_mtl_res, queries[query_id]);
    double duration_aqppp = (clock() - t_aqppp_start) / CLOCKS_PER_SEC;

    // Exact value... not necessary since I compute this elsewhere already?
    double t_exact_start = clock();
    double exact_value = expDemo::QueryRealValue(
        queries[query_id], PAR.TABLE_NAME, sql_connection_handle, PAR, "sum");
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
    time_sampling += duration_sampling;
    time_aqppp += duration_aqppp;
    time_exact += duration_exact;
    fprintf(results_file, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            query_id, result_sampling.first, result_sampling.second,
            result_aqppp.first, result_aqppp.second, exact_value, time_sampling,
            time_aqppp, time_exact, error_sampling, error_aqppp);

    n_queries_completed++;
  }
  fclose(results_file);

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

  fprintf(info_file, "Total experiment duration:   %f s\n", time_exp);
  fprintf(info_file, "Sample generation time:      %f s\n",
          time_samples_creation.first);
  fprintf(info_file, "Subsample generation time:   %f s\n",
          time_samples_creation.second);
  fprintf(info_file, "\n");

  fprintf(info_file, "Total sample rows:           %f\n",
          (double)sample[0].size());
  fprintf(info_file, "\n");

  fprintf(info_file, "Total pre-processing time:   %f s\n", time_preprocessing);
  fprintf(info_file, "Total pre-query time:        %f s\n",
          time_before_queries);
  fprintf(info_file, "Average query time sampling: %f s\n",
          time_sampling / n_queries_completed + time_read_sample);
  fprintf(info_file, "Average query time AQP++:    %f s\n",
          time_aqppp / n_queries_completed + time_read_sample +
              time_read_small_sample);
  fprintf(info_file, "Average query time exact:    %f s\n",
          time_exact / n_queries_completed);
  fprintf(info_file, "\n");

  fprintf(info_file, "Average error sampling:	  %f%%\n",
          100 * error_sampling_mean);
  fprintf(info_file, "Average error AQP++:		  %f%%\n",
          100 * error_aqppp_mean);
  fprintf(info_file, "Median error sampling:	  %f%%\n",
          100 * error_sampling_median);
  fprintf(info_file, "Median error AQP++:		  %f%%\n",
          100 * error_aqppp_median);
  fprintf(info_file, "99th percentile sampling:	  %f%%\n",
          100 * error_sampling_99p);
  fprintf(info_file, "99th percentile AQP++:	  %f%%\n",
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