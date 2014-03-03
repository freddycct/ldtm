#!/bin/bash

#cat test.txt | ./ldtm --num_topics=4 --prefix=output/test --iterations=2 --init_mu=0.5 --learning_rate=1.0 --record_logll=true

bzcat acmdl_sample_input.txt.bz2 | ./ldtm --num_topics=26 --prefix=output/acmdl --iterations=2 --init_mu=0.5 --learning_rate=10.0 --record_logll=true
./gc --time_series=output/acmdl_users.txt --interactions=acmdl_social_bef_input.txt --prefix=output/acmdl

