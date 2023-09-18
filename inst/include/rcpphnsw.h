//  RcppHNSW -- Rcpp bindings to hnswlib library for Approximate Nearest
//  Neighbors
//
//  Copyright (C) 2023 James Melville
//
//  This file is part of RcppHNSW
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RCPPHNSW_H
#define RCPPHNSW_H

#include <Rcpp.h>

#define HNSWLIB_ERR_OVERRIDE Rcpp::Rcerr

#include "hnswlib.h"

#endif // RCPPHNSW_H