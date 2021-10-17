#include "Clusterizer3x3.h"

#include <iostream>

void Clusterizer3x3::Init() { Reset(); }

void Clusterizer3x3::Reset() {

  for (int ic = 0; ic < mNcolumns; ic++) {
    for (int ir = 0; ir < mNrows; ir++) {
      setCell(ir, ic, 0, 0);
    }
  }

  for (int i = 0; i < nNclustersMax; i++) {
    addCluster(0., 0., 0, 0, 0, 0);
  }
  mNClusters = 0; // must be after the above loop
}

void Clusterizer3x3::Process() {

  FindCluster();

  if (mDebug > 0)
    std::cout << "Number of Clusters found : " << mNClusters << std::endl;

  if (mDebug > 1) {
    for (int iclus = 0; iclus < mNClusters; iclus++)
      std::cout << "clus#" << iclus << " : E=" << getClusterE(iclus)
                << " GeV, T=" << getClusterT(iclus)
                << " ns : row=" << getClusterRow(iclus)
                << ", col=" << getClusterCol(iclus)
                << ", mult=" << getClusterMult(iclus)
                << ", stCode= " << getClusterStatusCode(iclus) << std::endl;
  }
}

void Clusterizer3x3::FindCluster() {

  float Emax = -1;
  int colMax = -1;
  int rowMax = -1;

  for (int ir = 0; ir < mNrows; ir++) {
    if (ir < mRowRangeMin || ir > mRowRangeMax)
      continue;
    for (int ic = 0; ic < mNcolumns; ic++) {
      if (ic < mColRangeMin || ic > mColRangeMax)
        continue;

      if (isCellUsed(ir, ic) || getCellE(ir, ic) < mEleadingMin)
        continue;

      if (Emax < getCellE(ir, ic)) {
        Emax = getCellE(ir, ic);
        colMax = ic;
        rowMax = ir;
      }
    }
  }

  if (Emax < 0)
    return;

  float Eclus = 0.;
  int mult = 0;

  for (int ir = rowMax - 1; ir <= rowMax + 1; ir++) {
    if (ir < mRowRangeMin || ir > mRowRangeMax)
      continue;
    for (int ic = colMax - 1; ic <= colMax + 1; ic++) {
      if (ic < mColRangeMin || ic > mColRangeMax)
        continue;

      if (isCellUsed(ir, ic) || getCellE(ir, ic) < mEagg)
        continue;

      Eclus += getCellE(ir, ic);
      mult++;
      setCellAsUsed(ir, ic);
    }
  }

  int status = 0;

  addCluster(Eclus, getCellT(rowMax, colMax), rowMax, colMax, mult, status);

  FindCluster();
}
