#!/bin/bash


datum=05
day=14

Station=10410
#16546

./auswertung $(wc -l ${Station}.txt) ${datum} ${day} 00


