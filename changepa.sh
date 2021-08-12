#!/bin/bash

sed -i.bak "s/\(box.*\),[0-9.]*)/\1,${2}.000)/" regions/$1.reg
