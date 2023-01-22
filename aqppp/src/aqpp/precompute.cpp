// #include"stdafx.h"
#include "precompute.h"

#include <mutex>
#include <thread>

namespace aqppp {

Precompute::Precompute(std::string db_name, std::string table_name,
                       std::string agg_name,
                       std::vector<std::string> condition_names) {
  this->PAR.DB_NAME = db_name;
  this->PAR.AGGREGATE_NAME = agg_name;
  this->PAR.CONDITION_NAMES = condition_names;
  this->PAR.TABLE_NAME = table_name;
}

double Precompute::GetPrefixSumCube(
    const std::vector<std::vector<CA>>& mtl_points,
    SQLHANDLE& sql_connection_handle, MTL_STRU& o_mtl_res, std::string query_type,
    const DistId& distinct_ids) {
  MTL_STRU mtl_data = MTL_STRU();
  double t10 = clock();
  Precompute::ComputeDataCube(mtl_points, mtl_data, query_type, sql_connection_handle, distinct_ids);
  double mtl_data_time = (clock() - t10) / CLOCKS_PER_SEC;
  if (mtl_data.size() == 0) {
    o_mtl_res = MTL_STRU();
    return mtl_data_time;
  }
  double t1 = clock();
  const int DIM = mtl_points.size();
  std::vector<int> size_of_dimension = std::vector<int>(DIM);
  for (int i = 0; i < size_of_dimension.size(); i++) {
    size_of_dimension[i] = mtl_points[i].size();
  }
  std::vector<int> cur_ind = std::vector<int>(DIM);
  o_mtl_res = MTL_STRU(mtl_data);

  for (int i = DIM - 1; i >= 0; i--) {
    ComputeSumCube(0, i, size_of_dimension, cur_ind, o_mtl_res);
  }

  double run_time = (clock() - t1) / CLOCKS_PER_SEC + mtl_data_time;
  return run_time;
}

double Precompute::GetPrefixSumCube(const std::vector<int>& size_of_dimension,
                                    const MTL_STRU& mtl_data,
                                    MTL_STRU& o_mtl_res) {
  double start_time = clock();
  const int DIM = size_of_dimension.size();
  std::vector<int> cur_ind = std::vector<int>(DIM);
  o_mtl_res = MTL_STRU(mtl_data);
  for (int i = DIM - 1; i >= 0; i--) {
    ComputeSumCube(0, i, size_of_dimension, cur_ind, o_mtl_res);
  }
  double run_time = (clock() - start_time) / CLOCKS_PER_SEC;
  return run_time;
}

// cur_index need to be init outside with size of dimensions.
void Precompute::InitCube(int a, std::vector<int>& cur_indx,
                          MTL_STRU& o_mtl_data,
                          const std::vector<std::vector<CA>>& mtl_points) {
  if (a >= mtl_points.size()) {
    o_mtl_data.insert({cur_indx, 0});
    return;
  }
  for (int i = 0; i < mtl_points[a].size(); i++) {
    cur_indx[a] = i;
    InitCube(a + 1, cur_indx, o_mtl_data, mtl_points);
  }
  return;
}

int Precompute::MyLowerBound(const std::vector<CA>& cur_col, const CA& key) {
  int l = 0;
  int h = cur_col.size();  // Not n - 1
  while (l < h) {
    int mid = (l + h) / 2;
    if (DoubleLeq(key.condition_value, cur_col[mid].condition_value)) {
      h = mid;
    } else {
      l = mid + 1;
    }
  }
  return l;
}

void Precompute::ReadColumnOneThread(SQLHANDLE& sqlstatementhandle) {
  SQLFetch(sqlstatementhandle);
}

void Precompute::ComputeDataCubeOneThread(
    int thread_id, int start_row, int end_row, int n_conditions,
    const std::vector<std::vector<CA>>& mtl_points,
    double (&table)[BUFFER_COLULMS][BUFFER_ROWS], std::vector<double>& cube,
    std::vector<int>& cube_sizes, std::vector<std::mutex>& cube_locks) {
  std::vector<int> pos = std::vector<int>(n_conditions);
  // cout << GetCurrentProcessorNumber() << endl;
  // cout << thread_id << endl;
  for (int ri = start_row; ri <= end_row; ri++) {
    double acc_data = table[0][ri];
    for (int ci = 1; ci < n_conditions + 1; ci++) {
      CA cur_ca = CA();
      cur_ca.condition_value = table[ci][ri];

      int indi = MyLowerBound(mtl_points[ci - 1], cur_ca);
      pos[ci - 1] = indi;
    }
    int id = Tool::MapMultidimTo1Dim(pos, cube_sizes);
    if (id >= 0) {
      std::lock_guard<std::mutex> lock(cube_locks[id]);
      cube[id] += acc_data;
    }
  }
}

double Precompute::ComputeDataCube(
    const std::vector<std::vector<CA>>& mtl_points, MTL_STRU& o_mtl_data,
    const std::string query_type, const SQLHANDLE& sql_connection_handle, const DistId& distinct_ids) {
  // Setup
  double start_time = clock();
  const int THREAD_NUM_FOR_HANDLE = 3;
  int n_conditions = this->PAR.CONDITION_NAMES.size();
  assert(n_conditions < 10);   // else not enough memory allocated
  double table[BUFFER_COLULMS][BUFFER_ROWS] = {};  // DOES THIS NEED TO BE SO LARGE?
  o_mtl_data = MTL_STRU();
  
  // Compute cube size and initialise cube
  std::vector<int> cube_sizes = std::vector<int>();
  int total_cube_size = 1;
  for (int i = 0; i < mtl_points.size(); i++) {
    cube_sizes.push_back(mtl_points[i].size());
    total_cube_size *= mtl_points[i].size();
  }
  std::vector<double> cube = std::vector<double>(total_cube_size);
  std::vector<std::mutex> cube_locks = std::vector<std::mutex>(total_cube_size);

  // Define queries to be run
  std::vector<std::string> queries = std::vector<std::string>();
  queries.push_back("SELECT " + PAR.AGGREGATE_NAME + " FROM " + PAR.DB_NAME +
                    ".dbo." + PAR.TABLE_NAME + ";");
  for (auto cond_name : PAR.CONDITION_NAMES) {
    queries.push_back("SELECT CONVERT(varchar,  " + cond_name + ", 112) FROM " +
                      PAR.DB_NAME + ".dbo." + PAR.TABLE_NAME + ";");
  }  // TODO check this... 112 is for date columns... varchar probably doesn't make sense...




  // Initialise a SQL connection and statement handle for each query
  std::vector<SQLHANDLE> vec_sql_statement_handles =
      std::vector<SQLHANDLE>(n_conditions + 1);
   std::vector<SQLHANDLE> vec_sql_connections = std::vector<SQLHANDLE>();
    int t = 0;
  for (auto& statement_handle : vec_sql_statement_handles) {
     SQLHANDLE new_sql_connect_handle = NULL;
     SqlInterface::MakeSQLConnectionTrusted(PAR.DB_NAME, new_sql_connect_handle);
    // SqlInterface::MakeSQLConnectionWithPassword("D43139", "aqpplus", "aqpplus", new_sql_connect_handle);
     vec_sql_connections.push_back(new_sql_connect_handle);
     if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_STMT, new_sql_connect_handle,
                                       &statement_handle))
	//if (SQL_SUCCESS != SQLAllocHandle(SQL_HANDLE_STMT, sql_connection_handle,
 //                                     &statement_handle))
      return -1;
      t++;
  }





  std::vector<SQLULEN> pcrows = std::vector<SQLULEN>(n_conditions + 1);
  for (int ci = 0; ci < n_conditions + 1; ci++) {
    SqlInterface::SqlQuery(queries[ci], vec_sql_statement_handles[ci]);
    std::cout << "Query " << ci << " for computing data cube:\n"
              << queries[ci] << std::endl;
    SQLSetStmtAttr(vec_sql_statement_handles[ci], SQL_ROWSET_SIZE,
                   (void*)BUFFER_ROWS, 0);
    SQLBindCol(vec_sql_statement_handles[ci],  // Statement handle
               1,                              // Column number
               SQL_C_FLOAT,                   // C Data Type TODO need a single point of truth for this
               table[ci],                      // Data buffer
               BUFFER_ROWS,                  // Size of Data Buffer
               0);                             // Size of data returned
    SQLSetStmtAttr(vec_sql_statement_handles[ci], SQL_ATTR_ROW_ARRAY_SIZE,
                   (void*)BUFFER_ROWS, 0);
    SQLSetStmtAttr(vec_sql_statement_handles[ci], SQL_ATTR_ROWS_FETCHED_PTR,
                   &pcrows[ci], 0);
  }





  // cout << "cpu_num:"<<std::thread::hardware_concurrency() << endl;
  std::vector<std::thread> threads_for_read =
      std::vector<std::thread>(n_conditions);
  std::vector<std::thread> threads_for_handle =
      std::vector<std::thread>(THREAD_NUM_FOR_HANDLE);
  bool fetch_success = true;
  std::mutex fetch_lock;



  while (SQLFetch(vec_sql_statement_handles[0]) == SQL_SUCCESS) {
    for (int ci = 1; ci < n_conditions + 1; ci++) {
      threads_for_read[ci - 1] =
          std::thread(&Precompute::ReadColumnOneThread,
                      std::ref(vec_sql_statement_handles[ci]));
    }
    for (auto& th : threads_for_read) {
      th.join();
    }

    int start_row = 0;
    int end_row = -1;
    int left_rows = pcrows[0];
    for (int i = 0; i < THREAD_NUM_FOR_HANDLE; i++) {
      start_row = end_row + 1;
      end_row =
          fmin(pcrows[0] - 1, start_row - 1 +
                                  (left_rows + THREAD_NUM_FOR_HANDLE - i - 1) /
                                      (THREAD_NUM_FOR_HANDLE - i));
      if (end_row < start_row) break;
      left_rows -= (end_row - start_row + 1);
      threads_for_handle[i] = std::thread(
          &Precompute::ComputeDataCubeOneThread, i, start_row, end_row,
          n_conditions, std::ref(mtl_points), std::ref(table), std::ref(cube),
          std::ref(cube_sizes), std::ref(cube_locks));
    }
    for (auto& th : threads_for_handle) {
      th.join();
    }
    if (pcrows[0] < BUFFER_ROWS) break;
  }

  for (int id = 0; id < cube.size(); id++) {
    if (DoubleGreater(cube[id], 0)) {
      std::vector<int> indx = std::vector<int>();
      Tool::Map1DimToMultidim(id, cube_sizes, indx);
      o_mtl_data.insert({indx, cube[id]});
    }
  }

  double run_time = (clock() - start_time) / CLOCKS_PER_SEC;
  for (auto& handle : vec_sql_statement_handles) {
    SQLFreeHandle(SQL_HANDLE_STMT, handle);
  }
//   for (auto& connecthandle : vec_sql_connections) {
//     SQLDisconnect(connecthandle);
//     SQLFreeHandle(SQL_HANDLE_DBC, connecthandle);
//   }

  return run_time;
}

// cur_ind need to be init as a vector of len=dimensions outside
// o_mtl_res=mtl_data outside.
void Precompute::ComputeSumCube(int a, int b,
                                const std::vector<int>& size_of_dimension,
                                std::vector<int>& cur_ind,
                                MTL_STRU& o_mtl_res) {
  if (a == size_of_dimension.size()) {
    if (o_mtl_res.count(cur_ind) == 0) {
      o_mtl_res[cur_ind] = 0;
    }

    std::vector<int> pre_ind = std::vector<int>(cur_ind);
    pre_ind[b] -= 1;
    if (pre_ind[b] >= 0) o_mtl_res[cur_ind] += o_mtl_res[pre_ind];
    return;
  }

  for (int i = 0; i < size_of_dimension[a]; i++) {
    cur_ind[a] = i;
    ComputeSumCube(a + 1, b, size_of_dimension, cur_ind, o_mtl_res);
  }

  return;
}
}  // namespace aqppp