*** Settings ***
Documentation       Regression tests: run zCluster on a single cluster from 400SDSmall.csv (cluster 20,
...                 or cluster 1 for RubinDP2, as that is the only one in the DP2 footprint) for each
...                 supported database and compare against stored reference results.
...
...                 Run from the repository root with:
...                 robot --outputdir tests/results tests
...
...                 Useful options:
...                 --include SDSSDR12              run a single database (tags are database names)
...                 --exclude datalab              skip the (slow) NOIRLab Data Lab queries
...                 --exclude rubin                skip the Rubin queries (these are skipped anyway if there
...                                                is no RSP token in ~/.rsp-tap.token)
...                 --variable UPDATE_REFERENCE:True    re-write reference.json from this run
...                 --variable ZCLUSTER_ROOT:/path      test a different source tree
...
...                 Photometry is cached in tests/cache (via zCluster -c), so only the first run hits the network.
...                 Delete that directory to check that the remote services still work.
Library             ZClusterLibrary.py    ${ZCLUSTER_ROOT}
Test Template       Cluster Redshift Should Match Reference
Test Timeout        40 minutes


*** Variables ***
${ZCLUSTER_ROOT}        ${CURDIR}/..
${CATALOG}              ${ZCLUSTER_ROOT}/examples/400SDSmall.csv
${CACHE_DIR}            ${CURDIR}/cache
${REFERENCE_FILE}       ${CURDIR}/reference.json
${UPDATE_REFERENCE}     False
${Z_TOL}                0.02
${DELTA_REL_TOL}        0.1
${MAX_DZ_OVER_1PZ}      0.1


*** Test Cases ***    DATABASE    CLUSTER
SDSS DR12             SDSSDR12    20
    [Tags]    SDSSDR12

DECaLS DR10 (Data Lab)    DL_DECaLSDR10    20
    [Tags]    DL_DECaLSDR10    datalab

DECaLS DR11 (Data Lab)    DL_DECaLSDR11    20
    [Tags]    DL_DECaLSDR11    datalab

Rubin DP2             RubinDP2    1
    [Tags]    RubinDP2    rubin
    [Setup]    Skip Without RSP Token


*** Keywords ***
Cluster Redshift Should Match Reference
    [Arguments]    ${database}    ${cluster}
    ${fits}=    Run zCluster    ${CATALOG}    ${database}    ${cluster}
    ...    ${OUTPUT DIR}/runs/${database}    ${CACHE_DIR}/${database}
    ${result}=    Get Cluster Result    ${fits}    ${cluster}
    ${z_spec}=    Get Catalog Redshift    ${CATALOG}    ${cluster}
    Photo z Should Be Close To    ${result}    ${z_spec}    ${MAX_DZ_OVER_1PZ}
    Result Should Match Reference    ${result}    ${REFERENCE_FILE}    ${database}
    ...    ${Z_TOL}    ${DELTA_REL_TOL}    ${UPDATE_REFERENCE}

Skip Without RSP Token
    Skip If    not os.path.exists(os.path.expanduser('~/.rsp-tap.token'))
    ...    No RSP token in ~/.rsp-tap.token - see https://rsp.lsst.io/guides/auth/creating-user-tokens.html
