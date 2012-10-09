#!/bin/bash

[[ -e CorePyWrap.so ]] || lk CorePyWrap.so ~/ndb/cif-wrapper/wrapper/CorePyWrap.so
[[ -e bin/CorePyWrap.so ]] || lk bin/CorePyWrap.so ~/ndb/cif-wrapper/wrapper/CorePyWrap.so
