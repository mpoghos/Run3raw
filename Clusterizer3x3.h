#ifndef CLUSTERIZER3X3_H_
#define CLUSTERIZER3X3_H_

#include <Rtypes.h>
#include <tuple>
#include <vector>

const int mNcolumns = 48;
const int mNrows = 24;
const int nNclustersMax = 144; // 8*16 3x3 clusters possible

// Ecell, Tcell, IsUsed
typedef std::tuple<float, float, bool> Cell_t;
// Eclus, Tleading, RowL, ColL, Ncells, StatusCode
typedef std::tuple<float, float, int, int, int, int> Cluster_t;

class Clusterizer3x3 {

public:
  Clusterizer3x3() = default;
  ~Clusterizer3x3() = default;

  void Init();
  void Reset();

  void Process();
  void FindCluster();

  void setEaggMin(float e) { mEagg = e; }
  void setEleadingMin(float e) { mEleadingMin = e; }
  void setDebugLevel(int d) { mDebug = d; }

  void setRowRange(int imin, int imax) {
    mRowRangeMin = imin;
    mRowRangeMax = imax;
  }

  void setColRange(int imin, int imax) {
    mColRangeMin = imin;
    mColRangeMax = imax;
  }

  void setCell(int row, int col, Cell_t c) { mCells[row][col] = c; }
  void setCell(int row, int col, float acell, float tcell,
               bool status = false) {
    mCells[row][col] = std::make_tuple(acell, tcell, status);
  }

  void setCellAsUsed(int r, int c) { std::get<2>(mCells[r][c]) = true; }

  void addCluster(float e, float t, int row, int col, int mult, int status) {
    mClusters[mNClusters++] = std::make_tuple(e, t, row, col, mult, status);
  }

  Cell_t getCell(int row, int col) const { return mCells[row][col]; }
  float getCellE(int row, int col) const {
    return std::get<0>(mCells[row][col]);
  }
  float getCellT(int row, int col) const {
    return std::get<1>(mCells[row][col]);
  }
  bool getCellStatus(int row, int col) const {
    return std::get<2>(mCells[row][col]);
  }
  bool isCellUsed(int row, int col) { return getCellStatus(row, col); }

  int getNclusters() const { return mNClusters; }
  float getClusterE(int iclus) const { return std::get<0>(mClusters[iclus]); }
  float getClusterT(int iclus) const { return std::get<1>(mClusters[iclus]); }
  int getClusterRow(int iclus) const { return std::get<2>(mClusters[iclus]); }
  int getClusterCol(int iclus) const { return std::get<3>(mClusters[iclus]); }
  int getClusterMult(int iclus) const { return std::get<4>(mClusters[iclus]); }
  int getClusterStatusCode(int iclus) const {
    return std::get<5>(mClusters[iclus]);
  }

private:
  int mNClusters = 0;
  float mEagg = 0.1;
  float mEleadingMin = 0.15;
  int mDebug = -1;
  int mRowRangeMin = 0;
  int mRowRangeMax = mNrows;
  int mColRangeMin = 0;
  int mColRangeMax = mNcolumns;
  std::array<std::array<Cell_t, mNcolumns>, mNrows> mCells; //
  std::array<Cluster_t, nNclustersMax> mClusters;

  ClassDefNV(Clusterizer3x3, 1);
};

#endif