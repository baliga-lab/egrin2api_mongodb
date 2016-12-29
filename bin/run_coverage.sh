#!/bin/bash

PYTHONPATH=src EGRIN2API_SETTINGS=settings.cfg coverage run test/api_test.py && coverage xml --include=src/*
