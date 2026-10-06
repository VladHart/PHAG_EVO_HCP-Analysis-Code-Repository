# Analysis Code Repository and Controlled Data Access

## Overview

This repository contains the analysis code developed for the study. It is provided to support transparency, reproducibility, and review of the analytical methods.

The repository does **not** contain the dataset used in the analysis. No real, anonymized, pseudonymized, synthetic, sample, or demonstration dataset is provided through GitHub.

## Data availability

The data supporting the findings of this study are not publicly available in this GitHub repository due to applicable privacy, confidentiality, contractual, compliance, ethical, and regulatory requirements.

Access to the data may be considered following a reasonable and sufficiently detailed request to the corresponding author. Any request will be assessed in accordance with the relevant legal, ethical, contractual, compliance, data-protection, and regulatory frameworks, as well as any applicable approvals and restrictions governing the original data.

Subject to these requirements, access may require:

- a clearly defined scientific purpose and analysis plan;
- confirmation of the requester's identity, institutional affiliation, and relevant qualifications;
- approval by the data owner, study sponsor, or other responsible body, where applicable;
- execution of an appropriate data-use, confidentiality, or data-sharing agreement;
- evidence of ethics committee or institutional review board approval, where required;
- use of an approved secure analysis environment; and
- compliance with all applicable data-protection, information-security, contractual, and regulatory obligations.

Submission of a request does not guarantee that access will be granted. The scope and form of any access will depend on the applicable restrictions, permissions, and governance requirements. Where direct data sharing is not legally or contractually permitted, alternative approaches may be considered, where feasible and compliant.

### Data access requests

Reasonable inquiries should be addressed to:

**Corresponding author:** [Insert name]  
**Institution:** [Insert institution]  
**Email:** [Insert email address]  

Requests should include the proposed research question, scientific rationale, requested variables, intended analytical methods, planned outputs, requested access period, and information on data security and governance arrangements.

## Using the analysis code

The code contains one or more references to the location and/or name of the dataset used for the original analysis. These references are placeholders or environment-specific paths and will not provide access to the study data.

Before running the code, an authorized end user must replace each dataset reference with the actual local or secure-environment path and file name applicable to the authorized copy of the data.

For example, a placeholder reference such as:

```text
/path/to/authorized/dataset.ext
```

must be replaced on the end user's side with the actual approved location, for example:

```text
/approved/local/or/secure/path/actual_dataset.ext
```

The exact syntax will depend on the programming language, operating system, file format, and approved analysis environment. Users should search the code for dataset path variables, file-import commands, connection strings, or clearly marked placeholders such as:

```text
DATA_PATH
DATA_FILE
PATH_TO_DATA
REPLACE_WITH_AUTHORIZED_DATA_PATH
```

Only users who have obtained the necessary authorization and access permissions should configure the code to use the underlying data.

## No data included in this repository

To avoid ambiguity:

- no study dataset is stored in this repository;
- no dataset can be downloaded from this repository;
- no synthetic, sample, or demonstration dataset is provided;
- cloning or downloading this repository provides access only to the code and documentation; and
- the code cannot reproduce the study results unless an authorized user separately obtains the required data and correctly configures the dataset references.

## Responsibility of the end user

End users are responsible for ensuring that they:

- have all permissions required to access and analyze the data;
- use the data only for the approved purpose;
- work within an approved technical environment;
- do not upload confidential or protected data to GitHub;
- do not commit passwords, credentials, access tokens, connection strings, or other secrets to the repository; and
- comply with all applicable institutional, contractual, legal, ethical, compliance, information-security, and regulatory requirements.

## Reproducibility considerations

The code reflects the analytical workflow used for the study. Successful execution may depend on the authorized dataset having the required structure, variable names, coding conventions, formats, and data-processing history. Users may need to adapt file paths and environment-specific settings, but changes to the analytical logic should be documented to preserve transparency and traceability.

## Contact

For questions concerning the code, data-access process, or permitted use, please contact the corresponding author using the details provided above.
