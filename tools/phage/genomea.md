This report details issues with your submitted genome annotations. Your genome received a score of {{score}}

## Required Changes

The changes detailed in this section are required for acceptance of your submission.

### Missing Gene Features

These coding sequences are missing the associated gene feature. This is required for validation by NCBI's rules.

Feature ID | Location
---------- | ----------
{% for row in missing_genes %}
{{ row.id }} | `{{row.location.start}}..{{row.location.end}} [{{row.strand}}]` 
{% endfor %}

### Missing Product Tags

{{missing_tags_good}} / {{missing_tags_good + missing_tags_bad}}

Feature | Qualifiers
------- | ----------
{% for row in missing_tags %}
{{ row.id }} |
    {% for key in row.qualifiers %}
    {{ key }}
    <ul>
        {% for value in row.qualifiers[key] %}
        <li>{{value}}</li>
        {% endfor %}
    </ul>
    {% endfor %}
{% endfor %}


## Suggested Changes

These changes are not required, but are strongly encouraged in order to provide
a uniform genome annotation within the phage community.

### Start Codons

Start Codon | Count | Message
----------- | ----- | -------
{% for codon_key in weird_starts_overall_sorted_keys %}
{{ codon_key }} | {{ weird_starts_overall[codon_key] }} | {{ codon_key.__error }}
{% endfor %}

### Unannotated RBSs

Feature Type | ID | Location | Error | Upstream (-{{upstream_max}} .. -{{upstream_min}})
------------ | -- | -------- | ----- | -------------------------------------------------
{% for row in missing_rbs %}
{{ row.type }} | {{ row.id }} | `{{row.location.start}}..{{row.location.end}} [{{row.strand}}]` | {{row.__message}} | `{{row.__upstream}}`
{% endfor %}

## Notes

These areas may be indicative of a problem, or may simply be
informational. You should examine the areas mentioned in detail.

### Excessive Gaps

Region | Size | Bounding Gene Transcription Direction | Messages
------ | ---- | ------------------------------------- | --------
{% for row in excessive_gap %}
{{row[0]}} .. {{row[1]}} | {{row[1] - row[0]}} | {{row[2] | nice_strand}} {{row[3] | nice_strand}} |  {% if row[4] != 0 %}{{row[4]}} possible genes found in this region{% endif %}
{% endfor %}

### Excessive Overlaps

Feature A | Feature B | Shared Region | Overlap Length
--------- | --------- | ------------- | --------------
{% for row in excessive_overlap %}
{{row[0].id}} ({{row[0].location}}) | {{row[1].id}} ({{row[1].location}}) | {{row[2]}}..{{row[3]}} | {{row[3] - row[2]}}
{% endfor %}

### Coding Density

You have a coding density of {{ coding_density_real }}% which scores 
{{ coding_density }} / 100 on our scale. Most genomes should be in the 90% to 100%
coding density range
