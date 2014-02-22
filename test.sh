#!/bin/bash

cat test.txt | ./ldtm --num_topics=4 --prefix=output/test --iterations=2 --init_mu=0.5 --learning_rate=1.0 --record_logll=true
