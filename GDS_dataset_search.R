# Read the file line by line
txt_lines <- readLines("GEO_Datasets/Colorectal cancer.txt")

# Extract entries using regular expression
entries <- unlist(strsplit(paste(txt_lines, collapse = "\n"), "\\n\\s*\\d+\\.\\s+"))
entries <- entries[nchar(entries) > 0]  # remove empty entries

# Helper function to extract specific fields
extract_field <- function(entry, field) {
  pattern <- paste0(field, ":\\s*([^\n]+)")
  match <- regmatches(entry, regexpr(pattern, entry, perl = TRUE))
  if (length(match) == 0) return("")  # return empty if not found
  sub(paste0(field, ":\\s*"), "", match)
}


# Initialize empty list to hold results
parsed_data <- lapply(entries, function(entry) {
  list(
    Title = strsplit(entry, "\n")[[1]][1],
    Organism = extract_field(entry, "Organism"),
    Type = extract_field(entry, "Type"),
    Platform = extract_field(entry, "Platform[s]?"),
    Samples = sub(".*?([0-9]+) Samples.*", "\\1", entry),
    Accession = extract_field(entry, "Accession"),
    FTP = extract_field(entry, "FTP download")
  )
})

# Convert to data frame
df <- do.call(rbind, lapply(parsed_data, as.data.frame, stringsAsFactors = FALSE))

# Write CSV file to directory
write.csv(df, "GEO_Datasets/Colorectal cancer.csv", row.names = FALSE)
