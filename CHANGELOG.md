# Changelog

v1.0.0 - November 2026

# Added
- Function to round summary statisitcs in analysis table to 3 decimal places
- Addition of a unit test for this function
- CHANGELOG.md
- Docker files 


Unreleased - October 2026
# Added
- Added create_onyx_json to generate an Onyx analysis-fields JSON file from unclassified read-QC results.
- Added support for per-read summaries for both single-end and paired-end inputs.
- Added end-to-end tests for JSON generation with single-end and paired-end inputs.
- Added a Dockerfile for installing and running the Python application.
# Changed
- Updated the Onyx JSON-generation command to write the JSON file by default, rather than choosing between JSON storage and direct test/production Onyx writes.
- Improved paired-end read-count handling and the headline result.
- Return per-read results directly from the summary function.
- Updated package requirements and entry points for the Onyx JSON-generation tool.
- Expanded the README with command-line usage and example Onyx output.
- Refactored and cleaned up tests, including use of temporary test paths.
# Notes
- The FASTQ statistics tool’s version was updated to 0.2.0.
