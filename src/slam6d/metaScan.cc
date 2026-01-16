/*
 * metaScan implementation
 *
 * Copyright (C) by the 3DTK contributors
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief metascan is a collection of scans whcih can be treated just as
 *        a single scan
 *
 * @author Andreas Nuechter. Jacobs University Bremen gGmbH, Germany
 * @author Kai Lingemann. Inst. of CS, University of Osnabrueck, Germany.
 * @author Thomas Escher. Inst. of CS, University of Osnabrueck, Germany.
 */

#include "slam6d/metaScan.h"
#include "slam6d/kdMeta.h"

#ifdef WITH_METRICS
#include "slam6d/metrics.h"
#endif

MetaScan::MetaScan(std::vector<Scan*> scans, int nns_method) : m_scans(scans)
{
  // add this to the global vector for addFrame reasons
  Scan::allScans.push_back(this);

  if (m_scans.size() > 0) {
    m_identifier = "MetaScan("
      + std::string(m_scans.at(0)->getIdentifier())
      + "-" + std::string(m_scans.at(m_scans.size() - 1)->getIdentifier())
      + ")";
  } else {
    m_identifier = "metascan(empty)";
  }
}

MetaScan::~MetaScan()
{
  // remove this from the global vector
  for(ScanVector::iterator it = Scan::allScans.begin();
      it != Scan::allScans.end();
      ++it) {
    if(*it == this) {
      Scan::allScans.erase(it);
      break;
    }
  }
  // use clear to delete data
  for (auto& pair : m_pairs) {
    delete[] pair.second.first;
  }
  m_pairs.clear();
}

void MetaScan::createSearchTreePrivate()
{
#ifdef WITH_METRICS
  Timer tc = ClientMetric::create_metatree_time.start();
#endif //WITH_METRICS

  // TODO: there is no nns_type switch option for this one
  // because no reduced points are copied, this could be
  // implemented if e.g. cuda is required on metascans
  kd = new KDtreeMetaManaged(m_scans);

  /*
  Potential fix for above todo:
  Simliar to BasicScan, we now have m_pairs such that MetaScan
  is able to store its own points, if requested.
  So we might want to include some flag to indicate if we want to use
  a MetaManaged KDtree, or just use another type of tree, constructed from
  the MetaScan points stored in m_pairs (if available).
  */

  // DataXYZ xyz_orig(get("xyz reduced original"));
  // PointerArray<double> ar(xyz_orig);
  // switch(searchtree_nnstype)
  //   {
  //   case simpleKD:
  //     kd = new KDtree(ar.get(), xyz_orig.size(), searchtree_bucketsize);
  //     break;
  //   case ANNTree:
  //     kd = new ANNtree(ar, xyz_orig.size());
  //     break;
  //   case BOCTree:
  //     kd = new BOctTree<double>(ar.get(),
  //                               xyz_orig.size(),
  //                               10.0,
  //                               PointType(), true);
  //     break;
  //   case BruteForce:
  //       kd = new BruteForceNotATree(ar.get(),xyz_orig.size());
  //       break;
  //   case -1:
  //     throw std::runtime_error("Cannot create a SearchTree without setting a type.");
  //   default:
  //     throw std::runtime_error("SearchTree type not implemented");
  //   }

#ifdef WITH_METRICS
  ClientMetric::create_metatree_time.end(tc);
#endif //WITH_METRICS
}

size_t MetaScan::size() const
{
  return m_scans.size();
}

Scan* MetaScan::getScan(size_t i) const
{
  return m_scans.at(i);
}
