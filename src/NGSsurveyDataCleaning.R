## load libraries ----
library(tidyverse)
library(data.table)

## read data ----
# reading in the export of the google form for the survey
# the raw data is not available in the git - repository
data <- read_csv("data/Genomic Cancer NGS Report Elements 20221102.csv")


## clean data ----
# create shortcuts for column names
def.colnames <- tribble(
  ~colName, ~description,
  "timestamp", "Timestamp",
  "email", "Username",
  "profession", "What is your job title/position?",
  "region", "In which country do you work?",
  "organization", "What type of organization do you work for?",
  "cancer_type", "What kind of cancer cases  do you report/analyze?"
)

# rename columns 
setnames(data, def.colnames$description, def.colnames$colName)

# clean country names
data %>% 
  mutate(region = str_replace(region, "^U.*S.*", "USA"),
         region = str_replace(region, "^U.*K.*", "UK"),
         region = str_replace(region, "Brasil", "Brazil"), 
         region = str_replace(region, "Hong Kong", "China"), 
         region = str_replace(region, "China SAR", "China"),
         region = str_replace(region, "Espa.a", "Spain"),
         region = str_replace(region, "Deutschland", "Germany"))  -> data

## create IDs to distinguish between multiple responders ----
# remove emails for data security (not needed in this context)
data %>% 
  group_by(email) %>%
  mutate(ID = cur_group_id()) %>% 
  arrange(ID) %>% 
  ungroup() %>% 
  select(-email) -> data

## define groups for further analysis ----
# distinguish between Europe and North America
countries.europe <- c("Germany", "Belgium", "Spain", "UK", "Turkey")
countries.northamerica <- c("USA", "Canada")
data %>% 
  mutate(continent = ifelse(region %in% countries.europe, "Europe", "Other")) %>% 
  mutate(continent = ifelse(region %in% countries.northamerica, "North America", continent)) -> data

# create long-format
data %>% 
  pivot_longer(-c(timestamp, ID, region, continent)) -> data.long

# clean name column
data.long %>% 
  mutate(name = str_remove(name, "Which of these report elements are essential for a clinical report or not essential but preferred\\? \\(some of these may not be applicable to your practice\\) ")) -> data.long


# define groups of analysis
data.groups <- c("general", "molecular profile", "clinical trial")
data.long %>% 
  mutate(group = ifelse(grepl("molecular profile", name), "molecular profile", NA)) %>% 
  mutate(group = ifelse(grepl("[cC]linical [tT]rial", name), "clinical trial", group)) %>% 
  mutate(group = ifelse(grepl("^\\[.*\\]$", name), "general", group)) %>% 
  replace_na(list(group = "additional")) -> data.long

# define sections for the general information
data.long %>% 
  mutate(name = str_remove(name, "\\[")) %>% 
  mutate(name = str_remove(name, "\\]")) %>% 
  mutate(name = trimws(name)) %>% 
  mutate(name = str_remove(name, ".*\\?")) %>% 
  mutate(name = trimws(name)) -> data.long


names.overview <- c("Qualitative NGS Result Summary (i.e. This result is positive/negative)", "Diagnosis", 
                    "Tissue Source ", "Tumor purity/percentage", 
                    "Methods/Disclaimer: Assay limitations", "Methods/Disclaimer: Assay performance characteristics",
                    "Other disclaimers")
names.variants <- c("Reported list of genetic variants",
                    "Gene symbol in variant list",
                    "p. nomenclature in variant list",
                    "c. nomenclature in variant list",
                    "g. nomenclature in variant list",
                    "Gene transcript per variant",
                    "Variant allele frequency per variant")
names.functional <- c("Variant prevalence in patient's tumor type",
                      "Somatic variants classified by AMP/ASCO/CAP tiers",
                      "Somatic variants classified for pathogenicity/oncogenicity",
                      "Germline variants classified by ACMG/AMP criteria",
                      "Functional mechanism per variant (i.e. loss of function, activating)",
                      "Other Guidelines Used")
names.treatment <- c("Diagnostic information (if available) per variant",
                     "Prognostic information (if available) per variant",
                     "Therapeutic/Predictive information (if available) per variant",
                     "References used to interpret variant(s)",
                     "Databases/knowledgebases used to interpret variants",
                     "Versions of databases used to interpret variants",
                     "Clinical Trial Information")

data.long %>% 
  mutate(section = ifelse(name %in% names.overview, "Overview", NA)) %>% 
  mutate(section = ifelse(name %in% names.variants, "Variants", section)) %>% 
  mutate(section = ifelse(name %in% names.functional, "Functional", section)) %>% 
  mutate(section = ifelse(name %in% names.treatment, "Treatment", section))-> data.long

## separate responses for essentiality / report ----
responses.essential <- c("No response", "Not Recommended", "Preferred, non-essential", "Essential")
responses.reported <- c("Reported", "Not Reported", "No response")

data.long %>% 
  mutate(value = str_split(value, pattern = ";")) %>% 
  unnest(value) %>% 
  mutate(class = ifelse(value %in% responses.reported, "Report", NA)) %>% 
  mutate(class = ifelse(value %in% responses.essential, "Essentiality", class)) -> data.long

write_csv(data.long, "data/NGSreportElementsClean.csv")