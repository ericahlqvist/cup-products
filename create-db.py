import json
import psycopg2

with open('data/discriminants/p_rank_deg_[2_4_8]_4_mod_16.json', 'r') as file:
    data = json.load(file)

# Connect to the PostgreSQL database
conn = psycopg2.connect(
    database="cup-db",
    # user="Eric",
    # password="your_password",
    # host="your_host",
    # port="your_port"
)

# Create a database cursor
cursor = conn.cursor()
cursor_2 = conn.cursor()
# Iterate over your JSON data and insert into the table
for item in data:
    this_pol_txt = item['pol']
    #print(type(this_pol_txt))
    pol_txt = item['ext_pol'].strip('[|]')
    cyc_txt = item['ext_cyc'][1:-1]
    pol_vec = pol_txt.split(', ')
    cyc_vec = cyc_txt.split('], ')

    cursor_2.execute("""
            SELECT first_parent_pol, second_parent_pol FROM p_rank_deg_2_4_8 WHERE polynomial = %s;
        """, (this_pol_txt,))
    row = cursor_2.fetchone()
    print(row[0])
    print(row[1])
    for j in range(len(cyc_vec)-1):
        cyc_vec[j] = cyc_vec[j]+']'
    
    for i in range(len(pol_vec)):
        cursor.execute("""
            INSERT INTO p_rank_deg_2_4_16 (polynomial, first_parent_pol, second_parent_pol, third_parent_pol, cyc)
            VALUES (%s, %s, %s, %s, %s);
        """, (pol_vec[i], this_pol_txt, row[0], row[1], cyc_vec[i]))

# Commit the changes and close the connection
conn.commit()
cursor.close()
conn.close()
